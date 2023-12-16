# End-to-end Multi-sequence Alignment

# Install
1. Install Python 3.10;
2. Install dependencies via
```bash
pip3 install -r requirements.txt
```

# Run
1. Create a `./data` folder. Put all `.dna` files for all wild types into `./data/wt/`. Ensure that all file names are 
in the format of `<NAME>_sequence.dna`;
2. Put all sequencing results to be aligned into `./data/seq_results/`. Ensure that all files are in the format of 
`<NAME>-<GROUP_ID>-<ID>_PREMIX_<PLATE_NUM>.seq`;
3. (Optional) Put an Excel file named `Chromatogram_Report.xlsx` into `./data/seq_results`. This will allow the error 
message to show the trim value when the program failed to find the sub-sequence; 
4. Execute
```bash
python3 align.py
```

# Notes and Results
First, the program assumes that the primers stored in the `./data/wt/*.dna` are _in pairwise order_, i.e., all odd numbered 
primers are forward-facing, and all even numbered primers are in reverse.

The program will then attempt to find a sub-sequence that contains both primers. All such subsequences in the same group,
if exists, would be aligned using [muscle](https://github.com/rcedgar/muscle). The result will be stored as 
`<NAME>-<GROUP_ID>_<PRIMER_START>_<PRIMER_END>_alignment.html` in `./output`.