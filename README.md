# End-to-end Multi-sequence Alignment

# Install
1. For all commands, use Powershell for Windows and Terminal for macOS;
2. For macOS, install `brew` (a package manager that helps with installing and 
managing many modules) via
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
3. Install Git:
   - macOS: 
   ```bash
   brew install git
   ```
   - Windows: https://git-scm.com/download/win
   - Check your Git installation via
   ```bash
   git --version
   ```
4. Install Python 3.10:
   - macOS: https://www.python.org/downloads/release/python-31013/
   - Windows: https://www.microsoft.com/store/productId/9PJPW5LDXLZ5
   - Check your Python installation via
   ```bash
   python3 --version
   ```
   - Also make sure that Python and pip are in the same location, for macOS:
   ```bash
   which python3
   which pip3
   ```
   For Windows:
   ```bash
   gcm python3
   gcm pip3
   ```
   The two commands should return something in the same directory.
5. Download the script via
   ```bash
   cd <WHERE_YOU_WANT_TO_INSTALL>
   git clone https://github.com/fy-meng/seq-align.git
   cd seq-align
   ```
6. (Optional) Install `venv`. This is a virtual environment for Python. This 
allows the required packages to be installed in a standalone, virtual Python 
environment, and will not interfere with the dependencies of existing/future 
from other scripts. 
   1. Install `venv`:
      ```bash
      pip3 install venv
      ```
   2. Create an environment:
      ```bash
      python3 -m venv .venv
      ```
   3. Activate the environment. On macOS:
      ```bash
      source .venv/bin/activate
      ```
      On Windows:
      ```bash
      Set-ExecutionPolicy Unrestricted -Scope Process
      .venv\Scripts\activate
      ```
      You should now see a `(venv)` or `(.venv)` before the command line 
      prompt.
7. Install dependencies via
   ```bash
   pip3 install -r requirement.txt
   ```

# Usage
1. Create a `data` directory. Create two directories, `wt` and `seq_results` 
in `data`. You should have `./data/wt` and `./data/seq_results`. Put all `.dna` 
files into `./data/wt/`. All files names should be in the format of 
`<NAME>_sequence.dna`;
2. Put all `.seq` sequencing results to be aligned into `./data/seq_results/`. 
Ensure that all files are in the format of 
`<NAME>-<GROUP_ID>-<ID>_PREMIX_<PLATE_NUM>.seq`;
3. (Optional) Put an Excel file named `Chromatogram_Report.xlsx` into 
`./data/seq_results`. This will allow the error message to show the trim value 
when the program failed to find the sub-sequence; 
4. If `venv` is used during installation, activate it (if haven't already) in the 
same way as in [Install](#Install) 6. iii;
5. To run all DNAs in the `./data/wt/` directory:
```bash
python3 align.py
```
To run only selected DNAs:
```bash
python3 align.py <DNA1 DNA2 ...>
```

# Notes and Results
First, the program assumes that the primers stored in the `./data/wt/*.dna` are 
_in pairwise order_, i.e., all odd numbered primers are forward-facing, and all 
even numbered primers are in reverse.

The program will then attempt to find a sub-sequence that contains both 
primers. All such subsequences in the same group, if exists, would be aligned 
using [muscle](https://github.com/rcedgar/muscle). The result will be stored as 
`./output/<NAME>/<NAME>-<GROUP_ID>_<PRIMER_START>_<PRIMER_END>_alignment.html`.