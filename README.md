# End-to-end Multi-sequence Alignment

# Install

## Windows
For all commands, use Powershell.

1. Install Python 3.11 via Miniconda from https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html

   Follow the instructions from the graphical installer. All commands from now
   on is assumed to be run in Anaconda _Powershell_ Prompt.

   Follow the instructions from the graphical installer. Make sure that 
   `python` is linked to the Miniconda installation via:
   ```bash
   gcm python
   ```
   which should return `<MINICONDA_INSTALL>\python.exe`.

2. Install Git from https://git-scm.com/download/win

   Check your Git installation via
   ```bash
   git --version
   ```
3. Download the script via
   ```bash
   cd <WHERE_YOU_WANT_TO_INSTALL>
   git clone https://github.com/fy-meng/seq-align.git
   cd seq-align
   ```
   
## macOS
All commands are assumed to be run in Terminal.
1. Ensure that the default shell is zsh. Go to System Preferences - Users and 
   Groups - Control-click your name (not Command-click) - Advanced Options... -
   Login shell - change to `/bin/zsh`;
   
2. Install `brew` (a package manager that helps with installing and managing 
   modules) via
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
3. Install Git via
   ```bash
   brew install git
   ```
   Check your Git installation via
   ```bash
   git --version
   ```
4. Install Python via Miniconda from https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html

   __Important Note: Please install the Intel x86 version even if you are using an 
   Apple Silicon Mac.__

   Follow the instructions from the graphical installer. Make sure that 
   `python` is linked to the Miniconda installation via:
   ```bash
   which python
   ```
   which should return `<MINICONDA_INSTALL>/python`.
5. Install Muscle5 using conda via
   ```bash
   conda install bioconda::muscle
   ```
   Check your Muscle installation via
   ```bash
   muscle --version
   ```
6. Download the script via
   ```bash
   cd <WHERE_YOU_WANT_TO_INSTALL>
   git clone https://github.com/fy-meng/seq-align.git
   cd seq-align
   ```

# Usage
1. Open Anaconda Powershell Prompt (Windows) or Terminal (macOS). `cd` into the
   downloaded `seq-align` directory;
2. Activate the environment
   1. If this is the first time running the script, install dependencies. First, 
      create a conda virtual environment via
      ```bash
      conda create -n seq-align python=3.10
      ```
      Then activate the environment via
      ```bash
      conda activate seq-align
      ```
      You should see a `(seq-align)` before the prompt. Then install dependencies 
      via
      ```bash
      pip install -r requirements.txt
      ```
   2. Otherwise, activate the previously created environment via
   ```bash
   conda activate seq-align
   ```
   You should see a `(seq-align)` before the prompt.
3. Create a `data` directory. Create two directories, `wt` and `seq_results` 
in `data`. You should have `./data/wt` and `./data/seq_results`. Put all `.dna` 
files into `./data/wt/`. All files names should be in the format of 
`<NAME>_sequence.dna`;
4. Put all `.seq` sequencing results to be aligned into `./data/seq_results/`. 
Ensure that all files are in the format of 
`<NAME>-<GROUP_ID>-<ID>_PREMIX_<PLATE_NUM>.seq`;
5. (Optional) Put an Excel file named `Chromatogram_Report.xlsx` into 
`./data/seq_results`. This will allow the error message to show the trim value 
when the program failed to find the sub-sequence;
6. To run all DNAs in the `./data/wt/` directory:
```bash
python align.py
```
To run only selected DNAs:
```bash
python align.py <DNA1 DNA2 ...>
```

# Results
First, the program assumes that the primers stored in the `./data/wt/*.dna` are 
_in pairwise order_, i.e., all odd numbered primers are forward-facing, and all 
even numbered primers are in reverse.

The program will then attempt to find a sub-sequence that contains both 
primers. All such subsequences in the same group, if exists, would be aligned 
using [muscle](https://github.com/rcedgar/muscle). The result will be stored as 
`./output/<NAME>/<NAME>-<GROUP_ID>_<PRIMER_START>_<PRIMER_END>_alignment.html`.
