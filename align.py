import binascii
from collections import defaultdict
import copy
from math import ceil
import os
import platform
import re
import sys
from typing import *
import warnings

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import requests

WT_FILENAME_PATTERN = r'(.+?)_sequence.dna'
SEQ_FILENAME_PATTERN = r'(({tag}-\w+)[-\w+]+?)_.*\.seq'
PRIMER_PATTERN = r'<Primer recentID="\d+?" name="(.+?)" sequence="(\w+?)".*?><BindingSite location="(\d+)-(\d+)"'
PRIMER_SLICE_SIZE = 15
BASE_COLORS = {
    'A': (105, 248, 112),
    'T': (244, 159, 17),
    'C': (159, 205, 255),
    'G': (254, 230, 0)
}


def base2html(base: str) -> str:
    """
    Converts a single DNA base into HTML span with corresponding background
    color.
    :param base: the DNA base
    :return: an HTML span of the give base
    """
    if base in BASE_COLORS.keys():
        c = BASE_COLORS[base]
        return f'<span style="background-color: rgb({c[0]},{c[1]},{c[2]})">{base}</span>'
    else:
        return f'<span>{base}</span>'


def aln2html(records: List[SeqRecord], line_split=80):
    """
    Formats a list of aligned records into HTML.
    :param records: the list of aligned records.
    :param line_split: number of DNA bases per row.
    :return: the formatted HTML str
    """
    html = '<html><body><span style="font-family: Consolas,Menlo,monospace; font-size: 12px; display: inline-block;">\n'
    # handle ids
    ids = [r.id for r in records]
    max_id_len = max([len(s) for s in ids])
    max_id_len = 4 * (max_id_len // 4 + 1)
    ids = [s.ljust(4 * (max_id_len // 4 + 1)) for s in ids]

    # handle indices
    idx = [1 for _ in records]
    seq_len = len(records[0])
    max_idx_len = 4 * (len(str(seq_len)) // 4 + 1)

    for j in range(int(ceil(seq_len / line_split))):
        for i, r in enumerate(records):
            html += f'  <span style="white-space: pre;">{ids[i]}{str(idx[i]).ljust(max_idx_len)}</span>'
            for base in str(records[i].seq[j * line_split:(j + 1) * line_split]):
                html += base2html(base)
                if base != '-':
                    idx[i] += 1
            html += '  <br>\n'
        html += '  <br>\n'
    html += '</span></body></html>'
    return html


def findall(seq: Seq, pattern: Seq) -> List[int]:
    """
    Find all occurrences of a pattern in a DNA sequence.
    :param seq: the DNA sequence
    :param pattern: the pattern to be searched
    :return: list of indices of all occurrences
    """
    idx = []
    res = -1
    while True:
        res = seq.find(pattern, res + 1)
        if res == -1:
            break
        else:
            idx.append(res)
    return idx


def sub_seq(record: SeqRecord, start: str, end: str) -> Optional[SeqRecord]:
    """
    Slice a DNA sequence record between two patterns. Will attempt to search in
    both the forward and the reverse direction.
    :param record: the SeqRecord object containing the sequence
    :param start: the starting pattern
    :param end: the ending pattern
    :return: a SeqRecord object containing the subsequence if successful
    """
    start = Seq(start)
    end = Seq(end)

    # try forward search
    start_idx_f = findall(record.seq, start)
    end_idx_f = findall(record.seq, end)
    if len(start_idx_f) == 1 and len(end_idx_f) == 1:
        return record[start_idx_f[0]:end_idx_f[0] + len(end) + 1]

    # if failed, try search for reverse complement
    record = copy.deepcopy(record)
    record.seq = record.seq.reverse_complement()
    start_idx_r = findall(record.seq, start)
    end_idx_r = findall(record.seq, end)
    if len(start_idx_r) == 1 and len(end_idx_r) == 1:
        return record[start_idx_r[0]:end_idx_r[0] + len(end) + 1]

    # failed, raise warning
    err_msg = ''
    for direction, start_idx, end_idx in zip(['forward', 'reverse'],
                                             [start_idx_f, start_idx_r],
                                             [end_idx_f, end_idx_r]):
        err_msg += f'    {direction}: '
        if len(start_idx) == 0:
            err_msg += f'PRIMER_F not found; '
        elif len(start_idx) > 1:
            err_msg += f'PRIMER_F duplicated; '
        if len(end_idx) == 0:
            err_msg += f'PRIMER_R not found;'
        elif len(end_idx) > 1:
            err_msg += f'PRIMER_R duplicated;'
        err_msg += '\n'
    warnings.warn(err_msg)
    return None


def align_seqs(wt_filepath: str, seq_filepaths: List[str], primers: Iterable[Tuple[str, str, str, str]],
               primer_slice_size=PRIMER_SLICE_SIZE):
    """
    Align all sequences with the same tag, between PRIMER_START and PRIMER_END.
    The function will read the wild type files and the sequences and attempt to
    search for the primers. If successful, then align them using ClustalOmega,
    and then write the alignment result to an HTML file.
    :param wt_filepath: path to the .dna wild type file
    :param seq_filepaths: list of paths to the FASTA mutation files
    :param primers: list of (primer_start_name, primer_start_seq,
    primer_end_name, primer_end_seq) tuples
    :param primer_slice_size: number of bases to search for the primer. Only
    the innermost primer_slice_size number of bases will be searched in each
    sequence (the last few bases for the forward primers and the first few
    bases for the reverse primers).
    :param verbose: whether to raise a warning if the search failed.
    :return:
    """
    # read wild type files
    wt_filename = os.path.split(wt_filepath)[1]
    assert (m := re.match(WT_FILENAME_PATTERN, wt_filename))
    tag = m.group(1)

    with open(wt_filepath, 'rb') as handle:
        for wt in SeqIO.parse(handle, 'snapgene'):
            wt.id = f'{tag}_WT'

    # read mutation files
    seqs = defaultdict(lambda: dict())
    for seq_filepath in seq_filepaths:
        seq_filename = os.path.split(seq_filepath)[1]
        assert (m := re.match(SEQ_FILENAME_PATTERN.format(tag=tag), seq_filename))
        seq_name = m.group(2)

        with open(f'data/seq_results/{seq_filename}', 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seqs[seq_name][m.group(1)] = record

    for seq_name, d in seqs.items():
        success_dict = {seq_filename: False for seq_filename in d}
        err_msg_dict = {seq_filename: f'{seq_filename} failed:\n' for seq_filename in d}
        for primer_start_name, primer_start, primer_end_name, primer_end in primers:
            # move primers inwards
            primer_start = primer_start[-primer_slice_size:]
            primer_end = Seq(primer_end).reverse_complement()
            primer_end = primer_end[:primer_slice_size]

            # slice wild type
            wt_sub = sub_seq(wt, primer_start, primer_end)

            records = []
            for seq_filename, record in d.items():
                with warnings.catch_warnings(record=True) as w:
                    r = sub_seq(record, primer_start, primer_end)
                    if r:
                        records.append(r)
                if len(w) > 0 and issubclass(w[0].category, UserWarning):
                    err_msg = str(w[0].message)
                    err_msg = err_msg.replace('PRIMER_F', primer_start_name)
                    err_msg = err_msg.replace('PRIMER_R', primer_end_name)
                    err_msg_dict[seq_filename] += err_msg
                else:
                    success_dict[seq_filename] = True

            for seq_filename in d.keys():
                if not success_dict[seq_filename]:
                    warnings.warn(err_msg_dict[seq_filename])

            # output temporary multi-sequence file
            if len(records) > 0:
                if not os.path.isdir('tmp'):
                    os.mkdir('tmp')
                with open('tmp/input_clustal.fa', 'w+') as handle:
                    SeqIO.write(wt_sub, handle, 'fasta')
                with open('tmp/input_clustal.fa', 'a') as handle:
                    SeqIO.write(records, handle, 'fasta')

                # multi-sequence alignment using clustal-omega
                if not os.path.isdir('output'):
                    os.mkdir('output')
                if not os.path.isdir(os.path.join('output', tag)):
                    os.mkdir(os.path.join('output', tag))

                cmd = f'muscle -align tmp/input_clustal.fa -output tmp/alignment.aln -quiet'
                os.system(cmd)

                # convert to html
                with open('tmp/alignment.aln', 'r') as handle:
                    for alignment in AlignIO.parse(handle, 'fasta'):
                        alignment = sorted(alignment, key=lambda r: (not r.id.endswith('WT'), r.id))
                        html = aln2html(alignment)

                output_filepath = f'output/{tag}/{seq_name}_{primer_start_name}_{primer_end_name}_alignment.html'
                with open(output_filepath, 'w+') as handle:
                    handle.write(html)


def find_primers(file_name: str) -> List[Tuple[str, str]]:
    """
    Find all primers (names and DNA sequences) within a SnapGene DNA sequence.
    :param file_name: path to the SnapGene DNA file.
    :return: list of tuples of (primer_name, primer_seq, primer_start_idx, primer_end_idx)
    """
    with open(file_name, 'rb') as f:
        dna_data = binascii.hexlify(f.read())
        dna_data = bytes.fromhex(dna_data.decode())
        dna_data = str(dna_data)
    primers = list(re.findall(PRIMER_PATTERN, dna_data))
    primers = sorted(primers, key=lambda p: int(p[2]))  # sorted by their location
    return primers


def primer_combinations(primers_f: List[Tuple[str, str]], primers_r: List[Tuple[str, str]]) \
        -> List[Tuple[str, str, str, str]]:
    """
    Combine primers into pairs of forward and reverse primers. Assume that
    every other primer in the list has the same direction
    ([F1, R1, F2, R2, ...]), this function return all combinations of valid
    forward and reverse primer pairs
    ([(F1, R1), (F2, R2), ..., (F2, R2), (F2, R3), ...]).
    :param primers_f: list of tuples of (PRIMER_NAME, PRIMER_SEQ)
    :param primers_r: list of tuples of (PRIMER_NAME, PRIMER_SEQ)
    :return: list of tuples of (primer_start_name, primer_start_seq,
    primer_end_name, primer_end_seq)
    """
    result = []
    for f in primers_f:
        for r in primers_r:
            # if f.start_idx < r.start_idx
            if int(f[2]) < int(r[2]):
                result.append((f[0], f[1], r[0], r[1]))
    return result


def main():
    for wt_filename in os.listdir('./data/wt'):
        if m := re.match(WT_FILENAME_PATTERN, wt_filename):
            tag = m.group(1)
            # if no cmd args, run all
            # otherwise, only run the ones in the args
            if len(sys.argv) == 1 or tag in sys.argv:
                wt_filepath = os.path.join('./data/wt', wt_filename)

                # find primers
                primers = find_primers(os.path.join('./data/wt', wt_filename))
                primers_f = [p for p in primers if p[0].__contains__('F')]
                primers_r = [p for p in primers if p[0].__contains__('R')]
                primers = primer_combinations(primers_f, primers_r)

                # fetch the corresponding mutation files
                seq_filepaths = []
                for seq_filename in os.listdir('./data/seq_results'):
                    if re.match(SEQ_FILENAME_PATTERN.format(tag=tag), seq_filename):
                        seq_filepaths.append(os.path.join('./data/seq_results', seq_filename))

                print(f'aligning {tag}...', file=sys.stderr)
                align_seqs(wt_filepath, seq_filepaths, primers)


if __name__ == '__main__':
    # If Windows, download muscle executable
    if platform.system() == 'Windows' and not os.path.exists('./muscle.exe'):
        url = 'https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.win64.exe'
        print(f'downloading {url}...')
        r = requests.get(url, allow_redirects=True)
        with open('./muscle.exe', 'wb+') as f:
            f.write(r.content)

    # try open report excel
    try:
        report = pd.read_excel('./data/seq_results/Chromatogram_Report.xlsx', index_col=0)
    except FileNotFoundError:
        report = pd.DataFrame()

    # run
    main()
