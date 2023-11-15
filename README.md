# KEMET
**KE**gg **M**odule **E**valuation **T**ool  
_"KEMET - a python tool for KEGG Module evaluation and microbial genome annotation expansion"_  
[![DOI:10.1016/j.csbj.2022.03.015.svg](https://zenodo.org/badge/DOI/10.1016/j.csbj.2022.03.015.svg)](https://doi.org/10.1016/j.csbj.2022.03.015)

# Script description

The `kemet.py` script works as a command line tool that serves three main functions:

1) Evaluate **KEGG Modules Completeness** and summarize metabolic potential of MAGs/Genomes of interest, organizing the info into tables.  

2) Perform **HMM-based searches for ortholog genes (KO)** of interest, to expand KEGG Module Completeness evaluation.  

3) **Genome-scale models (GSMM or GEM) gapfill with evidence** from nucleotidic HMM searches, regarding KOs of interest.  

KEMET is well suited for metagenome-spanning analyses as well as single genomes usage, in order to get a better understanding of microbial metabolic and ecological functions.  

## Citing KEMET  
If you are using this program, please refer to our paper published in Computational and Structural Biotechnology Journal, available [here](https://doi.org/10.1016/j.csbj.2022.03.015).  


# Installation

The program is designed to have an easy installation procedure on UNIX-based machines, nonetheless the code is compatible with Windows systems.  

**Full installation is achieved in just a couple of minutes following the command-lines described in the [wiki pages](https://github.com/Matteopaluh/KEMET/wiki/0-Installation).**  

This tool was meant to have few external dependencies to ensure stability.  


## General Setup process and use - Conda environment

Refer to the [Setup wiki page](https://github.com/Matteopaluh/KEMET/wiki/1-Setup-process-using-a-Conda-environment) to properly set the working directory.  
Moreover it is important to follow the instructions to place relevant input files in the appropriate subdirectories and using proper format for said files.  

-----
# Command line (minimal required arguments)
```
./kemet.py [FASTA_file] -a [FORMAT] --hmm_mode [MODE] --gsmm_mode [MODE] (--skip_hmm) (--skip_gsmm) (--no_genome)
```

`[FASTA_file]`: FASTA file indication of the MAG/Genome of interest (with or without path indication e.g. `genomes/bin1.fasta`). With further arguments it can also be the indication of a KEGG annotation file.  

`-a [FORMAT]`: program used to annotate KEGG KOs, i.e. KEGG annotation format (either eggnog / kaas / kofamkoala) - used to generate KEGG MODULES recap tables. Default file extension must be maintained (e.g. `.emapper.annotations`, `.ko`)  

`--hmm_mode [MODE]`: when HMM analysis is desired, use this parameter to indicate a subset of KOs to search further using profile HMMs. `[MODE]` should be either one of `onebm`, `module`, `kos`, as described in the [wiki pages](https://github.com/Matteopaluh/KEMET/wiki).  

`--gsmm_mode [MODE]`: when GSMM/GEM gapfilling is desired, use this parameter to indicate whether to perform de-novo GSMM/GEM reconstruction or add reactions to an existing model. `[MODE]` should be either `denovo` or `existing`, respectively, as described in the [wiki pages](https://github.com/Matteopaluh/KEMET/wiki).  

`--skip_hmm`: use this to stop after KEGG MODULES Completeness Evaluation. The only output would be organized tables of metabolic potential.  

`--skip_gsmm`: use this to stop after HMM analysis.  

`--no_genome`: use this to indicate the (path to a) file with KEGG annotations, in order not to include MAG/Genome operations. Using this will result in stopping after KEGG MODULES Completeness Evaluation.  

-----
Other suggested optional parameters include:  

`--log`: store KEMET progress in a log file (STRONGLY suggested).  
`-v`: print more informations, for progress reporting purposes & more info.  
`-q`: print less informations and silence MAFTT and HMMER soft-errors (suggested).  

`--as_kegg`: changes how incomplete KEGG MODULES are summarized in recap tables - following KEGG-Mapper convention, i.e. Modules with less than 3 blocks are marked as INCOMPLETE regardless of the number of incomplete blocks. This imply using a more conservative approach regarding annotations.   

`--threshold_value [VALUE]`: use another quality filter to differentiate between legit HMM hits (default: 0.43).  

-----
## Script details  
For detailed info on the process/outputs of each KEMET task, as well as info on custom KEGG Modules & other, please refer to the [wiki pages](https://github.com/Matteopaluh/KEMET/wiki).  

-----
# Credits
Developed by Matteo Palù at Università degli Studi di Padova (2020-2023).