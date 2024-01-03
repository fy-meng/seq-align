import binascii
from collections import defaultdict
import copy
from math import ceil
import os
import platform
import re
import sys
import warnings

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import requests
from tqdm.auto import tqdm

PRIMER_PATTERN = r'<Primer recentID="\d+?" name="(.+?)" sequence="(\w+?)".*?><BindingSite location="(\d+)-(\d+)"'
PRIMER_SLICE_SIZE = 10
BASE_COLORS = {
    'A': (105, 248, 112),
    'T': (244, 159, 17),
    'C': (159, 205, 255),
    'G': (254, 230, 0)
}


def base2html(base):
    if base in BASE_COLORS.keys():
        c = BASE_COLORS[base]
        return f'<span style="background-color: rgb({c[0]},{c[1]},{c[2]})">{base}</span>'
    else:
        return f'<span>{base}</span>'


def aln2html(records, line_split=80):
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


def findall(seq: Seq, pattern: Seq):
    idx = []
    res = -1
    while True:
        res = seq.find(pattern, res + 1)
        if res == -1:
            break
        else:
            idx.append(res)
    return idx


def sub_seq(record: SeqRecord, start: str, end: str):
    start = Seq(start)
    end = Seq(end)

    # try forward search
    start_idx = findall(record.seq, start)
    end_idx = findall(record.seq, end)
    if len(start_idx) == 1 and len(end_idx) == 1:
        return record[start_idx[0]:end_idx[0] + len(end) + 1]

    # if failed, try search for reverse complement
    record = copy.deepcopy(record)
    record.seq = record.seq.reverse_complement()
    start_idx = findall(record.seq, start)
    end_idx = findall(record.seq, end)

    if len(start_idx) == 1 and len(end_idx) == 1:
        return record[start_idx[0]:end_idx[0] + len(end) + 1]

    # both failed
    return None


def align_seq(tag, primer_start, primer_end, primer_start_name, primer_end_name, verbose=False):
    # reading files
    with open(f'data/wt/{tag}_sequence.dna', 'rb') as handle:
        for wt in SeqIO.parse(handle, 'snapgene'):
            wt.id = f'{tag}_WT'

    # moving primers inwards
    primer_start = primer_start[-PRIMER_SLICE_SIZE:]
    primer_end = Seq(primer_end).reverse_complement()
    primer_end = primer_end[:PRIMER_SLICE_SIZE]

    mutations = defaultdict(lambda: dict())
    for file_name in os.listdir('data/seq_results'):
        if m := re.match(rf'(({tag}-\w+)[-\w]+?_PREMIX).*\.seq', file_name):
            seq_name = m.group(2)
            with open(f'data/seq_results/{file_name}', 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    mutations[seq_name][m.group(1)] = record

    # slicing
    wt = sub_seq(wt, primer_start, primer_end)
    for seq_name, d in mutations.items():
        records = []
        for file_name, record in d.items():
            r = sub_seq(record, primer_start, primer_end)
            if r:
                records.append(r)
            else:
                if verbose:
                    try:
                        trim = report["Trim"].loc[file_name]
                        warnings.warn(f'failed to find the subsequence for {record.name}, trim={trim}')
                    except KeyError:
                        warnings.warn(f'failed to find the subsequence for {record.name}')
        # output temporary multi-sequence file
        if len(records) > 0:
            if not os.path.isdir('tmp'):
                os.mkdir('tmp')
            with open('tmp/input_clustal.fa', 'w+') as handle:
                SeqIO.write(wt, handle, 'fasta')
            with open('tmp/input_clustal.fa', 'a') as handle:
                SeqIO.write(records, handle, 'fasta')

            # multi-sequence alignment using clustal-omega
            if not os.path.isdir('output'):
                os.mkdir('output')
            if not os.path.isdir(os.path.join('output', tag)):
                os.mkdir(os.path.join('output', tag))

            cmd = f'muscle -align tmp/input_clustal.fa -output tmp/alignment.aln -quiet'

            try:
                os.system(cmd)
            except Exception as e:
                warnings.warn(str(e))
                continue

            # convert to html
            with open('tmp/alignment.aln', 'r') as handle:
                for alignment in AlignIO.parse(handle, 'fasta'):
                    alignment = sorted(alignment, key=lambda r: (not r.id.endswith('WT'), r.id))
                    html = aln2html(alignment)

            with open(f'output/{tag}/{seq_name}_{primer_start_name}_{primer_end_name}_alignment.html', 'w+') as handle:
                handle.write(html)


def find_primers(file_name):
    with open(file_name, 'rb') as f:
        dna_data = binascii.hexlify(f.read())
        dna_data = bytes.fromhex(dna_data.decode())
        dna_data = str(dna_data)
    primers = list(re.findall(PRIMER_PATTERN, dna_data))
    primers = sorted(primers, key=lambda p: int(p[2]))
    assert len(primers) % 2 == 0
    primers = [primers[::2], primers[1::2]]
    result = []
    for i in range(len(primers[0])):
        for j in range(i, len(primers[1])):
            result.append({'primer_start_name': primers[0][i][0], 'primer_start': primers[0][i][1],
                           'primer_end_name': primers[1][j][0], 'primer_end': primers[1][j][1]})
    return result


def main():
    for file_name in os.listdir('./data/wt'):
        if m := re.match(r'(.+?)_sequence.dna', file_name):
            tag = m.group(1)
            # if no cmd args, run all
            # otherwise, only run the ones in the args
            if len(sys.argv) == 1 or tag in sys.argv:
                primers = find_primers(os.path.join('./data/wt', file_name))
                print(f'aligning {tag}...', file=sys.stderr)
                for p in tqdm(primers):
                    align_seq(tag, **p)


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
