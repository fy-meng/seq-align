from collections import OrderedDict
import os
import re
import traceback
from threading import Thread
from typing import *
import warnings

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter.scrolledtext import ScrolledText
from tkinter.ttk import Progressbar

from align import PRIMER_SLICE_SIZE, WT_FILENAME_PATTERN, SEQ_FILENAME_PATTERN
from align import find_primers, primer_combinations, align_seqs


def show_error(_self, *args):
    err = traceback.format_exception(*args)
    msg = ''.join(err)
    messagebox.showerror('Error', msg)


class MainView:
    MAX_STR_LEN = 60

    def __init__(self):
        self.root = tk.Tk()
        self.root.title('Multi-sequence Alignment Tool')
        self.root.minsize(width=600, height=180)

        self.primer_select_view: Optional[PrimerSelectView] = None

        self.header_size = tk.Label(self.root, text='Primer slicing size:')
        self.header_size.grid(row=0, column=0, padx=10, pady=5, sticky=tk.W)
        self.spinbox_size = tk.Spinbox(self.root, from_=0, to=30, width=3)
        self.spinbox_size.grid(row=0, column=2, padx=10, pady=5, sticky=tk.E)

        self.header_wt = tk.Label(self.root, text='Wild types:')
        self.header_wt.grid(row=1, column=0, padx=10, pady=5, sticky=tk.W)
        self.label_wt = tk.Label(self.root, text='', width=self.MAX_STR_LEN)
        self.label_wt.grid(row=1, column=1, padx=10, pady=5, sticky=tk.W)
        self.button_wt = tk.Button(self.root, text='Browse')
        self.button_wt.grid(row=1, column=2, padx=10, pady=5, sticky=tk.E)

        self.header_seq = tk.Label(self.root, text='Sequencing results:')
        self.header_seq.grid(row=2, column=0, padx=10, pady=5, sticky=tk.W)
        self.label_seq = tk.Label(self.root, text='', width=self.MAX_STR_LEN)
        self.label_seq.grid(row=2, column=1, padx=10, pady=5, sticky=tk.W)
        self.button_seq = tk.Button(self.root, text='Browse')
        self.button_seq.grid(row=2, column=2, padx=10, pady=5, sticky=tk.E)

        self.header_primers = tk.Label(self.root, text='Primers:')
        self.header_primers.grid(row=3, column=0, padx=10, pady=5, sticky=tk.W)
        self.label_primers = tk.Label(self.root, text='', width=self.MAX_STR_LEN)
        self.label_primers.grid(row=3, column=1, padx=10, pady=5, sticky=tk.W)
        self.button_primers = tk.Button(self.root, text='Select')
        self.button_primers['state'] = 'disabled'
        self.button_primers.grid(row=3, column=2, padx=10, pady=5, sticky=tk.E)

        self.button_submit = tk.Button(self.root, text='Submit')
        self.button_submit.grid(row=4, column=0, columnspan=3, padx=10, pady=5)

        self.progress_bar = None
        self.label_progress = None
        self.warning_msgbox = None

    def add_progress_info(self, progress_max_len: int):
        self.progress_bar = tk.ttk.Progressbar(self.root, maximum=progress_max_len)
        self.progress_bar.grid(row=5, column=0, columnspan=3, padx=10, pady=5, sticky='we')

        self.label_progress = tk.Label(self.root, height=2)
        self.label_progress.grid(row=6, column=0, columnspan=3, padx=10, pady=5, sticky=tk.W)

        self.warning_msgbox = ScrolledText(self.root)
        self.warning_msgbox.grid(row=7, column=0, columnspan=3, padx=10, pady=5, sticky=tk.W)

    def create_primer_select_window(self):
        self.primer_select_view = PrimerSelectView(self.root)

    def start_mainloop(self):
        self.root.mainloop()


class PrimerSelectView:
    def __init__(self, root):
        self.window = tk.Toplevel(root)
        self.window.title('Select Primers')
        self.window.minsize(width=450, height=400)

        self.frame_outer = tk.Frame(self.window)
        self.frame_outer.pack(padx=5, pady=5, fill=tk.BOTH, expand=True)
        self.frame_outer.columnconfigure(0, weight=1)
        self.frame_outer.rowconfigure(0, weight=1)

        self.canvas = tk.Canvas(self.frame_outer)
        self.canvas.grid(row=0, column=0, sticky='news')

        self.scrollbar = tk.Scrollbar(self.frame_outer, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scrollbar.grid(row=0, column=1, sticky='news')
        self.canvas.config(yscrollcommand=self.scrollbar.set)
        self.canvas.bind('<Configure>', lambda e: self.canvas.configure(scrollregion=self.canvas.bbox('all')))

        self.frame = tk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.frame, anchor='nw')

        self.headers = []
        self.labels = []
        self.listboxes = []
        self.listbox_scrollbars = []

        # in outer frame
        self.button = tk.Button(self.frame_outer, text='Select')
        self.button.grid(row=1, column=0, columnspan=2, pady=5)

    def add_list_boxes(self, wt_primers):
        for i, (wt_filename, primers) in enumerate(wt_primers.items()):
            header = tk.Label(self.frame, text=wt_filename)
            header.grid(row=3 * i, column=0, columnspan=4, padx=10, pady=5, sticky=tk.W)
            self.headers.append(header)

            _labels = []
            _listboxes = []
            _scrollbars = []

            for j, direction in enumerate(('Forward', 'Reverse')):
                # label
                label = tk.Label(self.frame, text=f'{direction} primers:')
                label.grid(row=3 * i + 1, column=2 * j, padx=10, pady=5, sticky=tk.W)
                _labels.append(label)

                # listbox
                listbox = tk.Listbox(self.frame, height=5, selectmode=tk.MULTIPLE, exportselection=0)
                listbox.grid(row=3 * i + 2, column=2 * j, padx=10, pady=5)
                # insert items and auto-select
                c = direction[0]  # 'F' or 'R'
                for k, (primer_name, _, _, _) in enumerate(primers):
                    listbox.insert(tk.END, primer_name)
                    if primer_name.__contains__(c):
                        listbox.select_set(k)
                _listboxes.append(listbox)

                # listbox scrollbar
                scrollbar = tk.Scrollbar(self.frame, orient=tk.VERTICAL, command=listbox.yview)
                scrollbar.grid(row=3 * i + 2, column=2 * j + 1, sticky='ns')
                listbox.config(yscrollcommand=scrollbar.set)
                scrollbar.config(command=listbox.yview)
                _scrollbars.append(scrollbar)

            self.labels.append(_labels)
            self.listboxes.append(_listboxes)
            self.listbox_scrollbars.append(_scrollbars)


class Controller:
    def __init__(self, view: MainView):
        self.view = view

        # variables
        self.primer_slicing_size = tk.IntVar(self.view.root, PRIMER_SLICE_SIZE)
        # stores both the wt_filenames and their primer info
        # keys are filenames of .dna wt files
        # items are list of (primer_name, primer_seq) tuples
        self.wt_primers_all: OrderedDict[str, List[Tuple[str, str]]] = OrderedDict()
        # keys are filenames of .dna wt files
        # items are list of (primer_start_name, primer_start, primer_end_name, primer_end) tuples
        self.wt_primers_pairs: OrderedDict[str, List[Tuple[str, str, str, str]]] = OrderedDict()
        self.seq_filepaths = []

        # bind variable
        self.view.spinbox_size.config(textvariable=self.primer_slicing_size)

        # bind functions
        self.view.button_wt.config(command=self.select_wt)
        self.view.button_seq.config(command=self.select_seq)
        self.view.button_primers.config(command=self.prompt_select_primers)
        self.view.button_submit.config(command=self.exec_threading)

    @staticmethod
    def format_filepaths(filepaths):
        # format paths to be relative, if possible
        filepaths = list(filepaths)
        cwd = os.path.abspath(os.getcwd())
        for i in range(len(filepaths)):
            filepaths[i] = os.path.abspath(filepaths[i])
            if filepaths[i].startswith(cwd):
                filepaths[i] = os.path.relpath(filepaths[i], start=os.getcwd())
        filepaths = sorted(filepaths)
        return filepaths

    def select_wt(self):
        filepaths = filedialog.askopenfilenames(initialdir=os.getcwd(),
                                                title='Select wild type files:',
                                                filetypes=[('SnapGene files', '*.dna'), ('All files', '*.*')])
        filepaths = self.format_filepaths(filepaths)

        # find all primers
        for filepath in filepaths:
            primers = find_primers(filepath)
            self.wt_primers_all[filepath] = primers

        # update label
        label_str = ', '.join(filepaths)
        if len(label_str) > self.view.MAX_STR_LEN:
            label_str = label_str[:self.view.MAX_STR_LEN - 3] + '...'
        self.view.label_wt['text'] = label_str

        # if success, enable the button
        if len(filepaths) > 0:
            self.view.button_primers['state'] = 'normal'
        else:
            self.view.button_primers['state'] = 'disabled'

        # then choose default primers
        label_strs = []
        for wt_filename, primers in self.wt_primers_all.items():
            primers_f = [p for p in primers if p[0].__contains__('F')]
            primers_r = [p for p in primers if p[0].__contains__('R')]
            self.wt_primers_pairs[wt_filename] = primer_combinations(primers_f, primers_r)
            primer_names_f = [p[0] for p in primers_f]
            primer_names_r = [p[0] for p in primers_r]
            label_strs.append(f'{wt_filename}: [({", ".join(primer_names_f)}), ({", ".join(primer_names_r)})]')
        label_str = ', '.join(label_strs)
        if len(label_str) > self.view.MAX_STR_LEN:
            label_str = label_str[:self.view.MAX_STR_LEN - 3] + '...'
        self.view.label_primers['text'] = label_str

    def select_seq(self):
        filepaths = filedialog.askopenfilenames(initialdir=os.getcwd(),
                                                title='Select sequencing results:',
                                                filetypes=[('FASTA files', '*.seq *.fasta *.fas *.fa'),
                                                           ('All files', '*.*')])
        filepaths = self.format_filepaths(filepaths)
        self.seq_filepaths.extend(filepaths)

        # update label
        label_str = ', '.join(filepaths)
        if len(label_str) > self.view.MAX_STR_LEN:
            label_str = label_str[:self.view.MAX_STR_LEN - 3] + '...'
        self.view.label_seq['text'] = label_str

    def prompt_select_primers(self):
        self.view.create_primer_select_window()
        self.view.primer_select_view.add_list_boxes(self.wt_primers_all)
        self.view.primer_select_view.button.config(command=self.select_primers)

    def select_primers(self):
        listboxes: List[Tuple[tk.Listbox, tk.Listbox]] = self.view.primer_select_view.listboxes

        label_strs = []

        for i, (wt_filename, (listbox_f, listbox_r)) in enumerate(zip(self.wt_primers_all.keys(), listboxes)):
            indices_f = listbox_f.curselection()
            primers_f = [self.wt_primers_all[wt_filename][j] for j in indices_f]
            indices_r = listbox_r.curselection()
            primers_r = [self.wt_primers_all[wt_filename][j] for j in indices_r]
            self.wt_primers_pairs[wt_filename] = primer_combinations(primers_f, primers_r)

            primer_names_f = [p[0] for p in primers_f]
            primer_names_r = [p[0] for p in primers_r]
            label_strs.append(f'{wt_filename}: [({", ".join(primer_names_f)}), ({", ".join(primer_names_r)})]')

        label_str = ', '.join(label_strs)
        if len(label_str) > self.view.MAX_STR_LEN:
            label_str = label_str[:self.view.MAX_STR_LEN - 3] + '...'
        self.view.label_primers['text'] = label_str

        self.view.primer_select_view.window.destroy()
        del self.view.primer_select_view

    def exec_threading(self):
        self.view.add_progress_info(len(self.wt_primers_pairs))
        t1 = Thread(target=self.exec)
        t1.start()

    def exec(self):
        for i, wt_filepath in enumerate(self.wt_primers_pairs.keys()):
            self.view.progress_bar.step()

            wt_filename = os.path.split(wt_filepath)[1]
            if m := re.match(WT_FILENAME_PATTERN, wt_filename):
                tag = m.group(1)
                self.view.label_progress['text'] = f'aligning {tag}...'

                _seq_filepaths = []
                for seq_filepath in self.seq_filepaths:
                    seq_filename = os.path.split(seq_filepath)[1]
                    if re.match(SEQ_FILENAME_PATTERN.format(tag=tag), seq_filename):
                        _seq_filepaths.append(os.path.join('./data/seq_results', seq_filepath))

                with warnings.catch_warnings(record=True) as ws:
                    align_seqs(wt_filepath, _seq_filepaths, self.wt_primers_pairs[wt_filepath])
                    for w in ws:
                        if issubclass(w.category, UserWarning):
                            self.view.warning_msgbox.insert(tk.INSERT, str(w.message))
                            self.view.warning_msgbox.see(tk.END)

        self.view.progress_bar.configure(mode="determinate", value=len(self.wt_primers_pairs))
        self.view.label_progress['text'] = 'finished'

    def start(self):
        self.view.start_mainloop()


if __name__ == '__main__':
    tk.Tk.report_callback_exception = show_error

    view = MainView()
    controller = Controller(view)
    controller.start()
