
# [<samp>Jousset2026 structures</samp>](https://github.com/desiro/Jousset2026_structures)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-bd0000.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python v3.9.7](https://img.shields.io/badge/Language-Python_v3-75a8d3.svg)](https://www.python.org/)
[![Conda v4.11.0](https://img.shields.io/badge/Uses-Conda-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Description

RNA structure prediction script for <samp>Jousset2026</samp>.

### Mandatory Prerequisites

* [![Python v3.9.7](https://img.shields.io/badge/Python_v3.9.7-75a8d3.svg)](https://www.python.org/downloads/release/python-397/)
* [![ViennaRNA v2.5.0](https://img.shields.io/badge/ViennaRNA_v2.5.0-006795.svg)](https://www.tbi.univie.ac.at/RNA/)

### Optional Prerequisites

* [![Conda v4.11.0](https://img.shields.io/badge/Conda_v4.11.0-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Installation

I recommend using Miniconda and following the steps below. If this is the first time using conda, you should probably restart your shell after the installation of Miniconda. The following will demonstrate the installation and set up of Miniconda on Linux, which should be similar on other platforms. For Windows 10 users, I advise using the [Ubuntu 20.04 LTS](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71?cid=msft_web_chart) subsystem. More information can be found on the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [Bioconda](https://bioconda.github.io/) pages.

### Conda Installation

Installing Miniconda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Updating Miniconda and setting channels:
```
conda update conda
conda update python
conda config --add channels bioconda
conda config --add channels conda-forge
```

Installing Conda packages:
```
git clone https://github.com/desiro/Jousset2026_structures.git
cd Jousset2026_structures
conda env create -f envStructuresSHAPE.yml
```

***

## Execution

```
conda run -n structures_SHAPE python getStructuresSHAPE.py --pickleData \
    --prefix <prefix> \
    --genomeFasta <genome.fa> --dataSHAPE <shape dir>

# <prefix>    -> output directory and prefix for result files
# <genome.fa> -> input genome FASTA file
# <shape dir> -> directory with all SHAPE-MaP files; names have to match the fasta headers
```

```
conda run -n structures_SHAPE python compareStructuresSHAPE.py --pickleData \
    --prefix <prefix> \
    --genomeFasta1 <genome 1.fa> --dataSHAPE1 <shape 1 dir> --structureFile1 <structures 1.pcl> \
    --genomeFasta2 <genome 2.fa> --dataSHAPE2 <shape 2 dir> --structureFile2 <structures 2.pcl>

# <prefix>           -> output directory and prefix for result files
# <genome 1.fa>      -> first input genome FASTA file
# <genome 2.fa>      -> second input genome 2 FASTA file
# <shape 1 dir>      -> first directory with all SHAPE-MaP files; names have to match the fasta headers
# <shape 2 dir>      -> second directory with all SHAPE-MaP files; names have to match the fasta headers
# <structures 1.pcl> -> first pickled file from getStructuresSHAPE
# <structures 2.pcl> -> second pickled file from getStructuresSHAPE
```

***

## Results

```
bash structuresSHAPE.sh
```

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite <samp>Jousset2026</samp>.

```
A.-C. Jousset, A. Hache, C. Jakob, B. Chane-Woon-Ming, D. Ferhadian, D. Desirò, R.P. Stansilaus, M. Marz, M. Schwemmle, H. Bolte and R. Marquet.
"The wild-type packaging signal network cooperatively induces a metastable conformation of the influenza A virus genome during packaging."
In review, 2026.
```
