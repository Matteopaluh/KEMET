# KEMET
**KE**gg **M**odule **E**valuation **T**ool  
_"KEMET - a python tool for KEGG Module evaluation and microbial genome annotation expansion"_  
[![DOI:10.1016/j.csbj.2022.03.015.svg](https://zenodo.org/badge/DOI/10.1016/j.csbj.2022.03.015.svg)](https://doi.org/10.1016/j.csbj.2022.03.015.svg)

## Script description

The `kemet.py` script works as a command line tool that serves three main functions:

1) Evaluate KEGG Modules Completeness and summarize the metabolic potential of MAGs/Genomes of interest, generating organized tables.  
2) Perform HMM-based searches for ortholog genes (KO) of interest, to expand the KEGG Module Completeness evaluation.
3) Genome-scale model gapfill with evidence from nucleotidic HMM searches, regarding the KOs of interest.

KEMET is well suited for single genomes usage as well as metagenome-spanning analyses, in order to get a better understanding of microbial metabolic and ecological functions. 

## Citing KEMET  
If you are using this program, please refer to our paper published in Computational and Structural Biotechnology Journal, available [here](https://doi.org/10.1016/j.csbj.2022.03.015).  

-----
## Installation

The program is designed to have an easier installation procedure on UNIX-based machines, nonetheless the code is compatible with Windows systems as well.  
(tested on: Ubuntu 20.04 LTS with Linux 5.13 kernel - Windows 10 build 19043.1526 - March 2022)  
Windows systems can use the "Windows Subsystem for Linux" extra feature, as described [here](https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html) by MAFFT developers.  

The easiest way to install KEMET using both UNIX or Windows subsystem is through the Anaconda package manager. You could see Anaconda documentation at [this link](https://docs.anaconda.com/).  

For Windows subsystems, a quick installation is achieved with:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

1) Create a new environment and install dependencies

With Anaconda properly set, it is possible to create an environment called "kemet" using
```
conda create --name kemet python=3.7
conda activate kemet
```  
To install other dependencies use the following commands:
```
conda install -c bioconda mafft
conda install -c bioconda hmmer
pip install reframed

# OPTIONAL (only if GSMM mode use is intended)

pip install carveme==1.4.1
conda install -c bioconda diamond
conda install -c bioconda prodigal
```

2) Clone the git repository in your preferred working directory, and access it.  
It's also possible to download a compressed version directly from this GitHub repository.

```
git clone https://github.com/Matteopaluh/KEMET.git
cd KEMET
```
-----
## General Use - Conda environment

0) Open a shell in the KEMET directory. Scripts execution should be enabled. If the opposite is true use  
`chmod +x ./*.py`

1) Run once the `setup.py` script (if genome-scale models functionalities are wanted, add the `-G` parameter).  
This will set the folders in which different input and outputs will be stored.

2) Set input files into proper paths (**IMPORTANT**):  
- Copy MAG/Genome sequences to be analysed in `KEMET/genomes/` folder, which is created after the setup process.  
**NOTE:**  
**Sequence file extensions are supposed to be ".fa",".fna" or ".fasta"**  
**FASTA header are supposed not to be repeated in a single MAG/Genome**.  
If necessary, rename MAGs/Genomes and FASTA headers accordingly.  

- Copy KEGG KOs annotations (derived from different sources) in `KEMET/KEGG_annotations/` folder, created in the setup process.  
The script requires an indication of the program used to generate input KEGG annotation (eggNOG, KofamKOALA, and KAAS or KAAS-like format are supported up to March 2022).  
**Annotations format don't need to be changed from their original output** (truncated example files can be found in `KEMET/toy/` folder)  

- If genome-scale model expansion is needed, optional pre-existing genome-scale models within the BiGG namespace (".xml" files) can be copied in the `KEMET/models/` folder, created after the setup process.  

## IMPORTANT NOTES:  
The **files' extension does not need to be modified**, neither the rest of file name, except **if** it does not refer to the same genome in input:  

e.g. **bin1.fasta/.fa/.fna** needs to be selected to check annotations from file **bin1.emapper.annotations**, and work with a genome-scale model **bin1.xml**.  

`KEMET/KEGG_MODULES/` folder presence as in GitHub is **necessary** for script usage. It contains files that represent KEGG Modules block structure (REF: [KEGG MODULE resource](https://www.genome.jp/kegg/module.html)); missing KO orthologs are deduced from these structures. Other "custom" Modules could be added to that folder, if formatted in the proper way.


3) (**ONLY if HMM and GSMM steps are needed**) Fill in the textual intruction file "genomes.instruction":  
excluding the header, each line should have a **tab-separated** indication of:

|MAG/Genome FASTA|Taxonomic indication|Metabolic model universe|
|---|---|---|  

- The MAG/Genome FASTA indicate the MAG/Genome of interest file name (e.g. bin1.fasta)  

- The taxonomic indication should be taken from the KEGG Brite taxonomic indication (specifically from the C-level, that most of the times coincide with NCBI phylum level taxonomy) (REF: [BRITE Organism table](https://www.kegg.jp/brite/br08601)) (e.g. Actinobacteria)  

- Metabolic model universe comprehend grampos, gramneg, archaea or other custom universe (this is an optional indication needed for GSMM de-novo reconstruction)  

4) (**ONLY if HMM and GSMM steps are needed**) Fill in other instruction text files.  
If HMM-analyses are desired, these need either the "module_file.instruction" or the "ko_file.instruction" files as follows, depending on the desired MODE OF USE (which needs to be specified with the `--hmm_mode MODE` parameter).  

|MODE|Analysis|Instructions|
|---|---|---|
|onebm|KOs from KEGG Modules missing 1 block|(No need to fill instruction files)|
|modules|KOs from a fixed list of KEGG Modules|(One per line indication in the "module_file.instruction" file)|
|kos|KOs from a fixed list of orthologs|(One per line indication in the "ko_file.instruction" file)|  

5) Launch the `kemet.py` command line script with your arguments of choice!

-----
# Command line (minimal, only required arguments)
```
./kemet.py [FASTA_file] -a [FORMAT] (--skip_hmm) --hmm_mode [MODE] (--skip_gsmm) --gsmm_mode [MODE]
```

`[FASTA_file]`: indication of the MAG/Genome of interest FASTA (with or without path indication e.g. `genomes/bin1.fasta`)  

`-a [FORMAT]`: indication of which program was used to annotate KEGG KOs, i.e. the KEGG annotation format (either eggnog/kaas/kofamkoala) to generate KEGG MODULES recap tables. Default file extension must be maintained (".emapper.annotations", ".ko")  

`--skip_hmm`: use this parameter to stop after KEGG MODULES Completeness Evaluation. The only output would be organized tables of metabolic potential obtained from previous annotation softwares.  

`--hmm_mode [MODE]`: if the HMM analysis is desired, this parameter is required to indicate which subset of KOs to search with profile HMMs (`onebm, module, kos`), which are further described.  

`--skip_gsmm`: use this parameter to stop after HMM analysis.  

`--gsmm_mode [MODE]`: if GSMM operations with HMM-hits are desired, this parameter is required to indicate whether to perform de-novo GSMM reconstruction or gapfilling, with `[MODE]` equal to `denovo` or `existing`, respectively, as further described.  

-----
Other optional parameters are described in the following paragraph:  

`-v`: print more informations, for progress reporting purposes & more info.  
`-q`: print less informations and silence MAFTT and HMMER soft-errors (suggested).  
`--log`: store KEMET progress in a log file (suggested).  

`--as_kegg`: change the way incomplete KEGG MODULES are considered in the metabolic potential recap tables, following KEGG-Mapper convention, i.e. Modules with less than 3 blocks that have at least 1 incomplete block are marked as INCOMPLETE regardless.  

`--threshold_value [VALUE]`: use a different quality filter to differentiate between HMM hits (default: 0.43).  

`--update_taxonomy_codes`: when downloading a new BRITE taxonomy with `setup.py`, this parameter is necessary to have an updated species (BRITE E-class) list for each C-class, that is to have an updated scope for taxonomic-informed HMM creation.  

`--skip_nt_download`: avoid the KEGG GENES sequences download, e.g. if they were previously downloaded, in order not to follow KEGG updates.  
`--skip_msa_and_hmmbuild`: avoid MAFFT and HMMER hmmbuild commands launch, e.g. if multi-alignment and pHMMs were done generated.  
`--retry_nhmmer`: use this parameter after a full pipeline run, to re-run nHMMER commands properly, e.g. in order to use a different scoring threshold.  

-----
# Script details  
For detailed info on the process or the outputs of each KEMET task, please refer to the [wiki page](https://github.com/Matteopaluh/KEMET/wiki)

-----
# Credits
Developed by Matteo Palù at Università degli Studi di Padova (2020-2022).