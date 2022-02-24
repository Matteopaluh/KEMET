# KEMET
**KE**gg **M**odule **E**valuation **T**ool  

## Script description

The `kemet.py` script works as a command line tool that serves three main functions:

1) Evaluate KEGG Modules Completeness, producing organized tables that summarize metabolic potential of MAGs/Genomes of interest. 
2) Perform HMM-based searches for ortholog genes (KO) of interest chosen from a subset after KEGG Module Completeness evaluation.
3) Genome-scale model gapfill with evidence from nucleotidic HMM searches, regarding KOs of interest.

It is well suited for single genomes usage as well as metagenome-spanning analyses, in order to get a better understanding of microbial metabolic and ecological functions.  

_"KEMET - a python tool for KEGG Module evaluation and microbial genome annotation expansion"_ - Manuscript under revision  

-----
## Installation

The program is designed to have an easier installation procedure on UNIX-based machines, nonetheless the code is compatible with Windows systems as well (tested on: Ubuntu 19.04 LTS - Windows 10 build 19043.1526 - February 2022).  
Windows systems can use extra features such as the "Windows Subsystem for Linux", as described [here](https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html) by MAFFT developers.  

The easiest way to install KEMET using both UNIX or Windows subsystem is to use the Anaconda package manager. You could see Anaconda documentation at [this link](https://docs.anaconda.com/).  

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

0) Open a shell in the KEMET directory. Scripts execution should be enabled. If the opposite is true use `chmod +x ./*.py`

1) Run once the `setup.py` script (if genome-scale models possibilities are wanted, add the `-G` parameter).  
This will set the folders in which different input and outputs will be stored.

2) Set input files into proper paths (**IMPORTANT**):  
- Copy MAG/Genome sequences to be analysed in `KEMET/genomes/` folder, which is created after the setup process.  
**NOTE:**  
**File extensions are supposed to be ".fa",".fna" or ".fasta"**  
**FASTA header are supposed not to be repeated in a single MAG/Genome**.  
If necessary, rename MAGs/Genomes and FASTA headers accordingly.  

- Copy KEGG KOs annotations (derived from different sources) in `KEMET/KEGG_annotations/` folder, created in the setup process.  
The script requires an indication of the program used to generate input KEGG annotation (eggNOG, KofamKOALA, and KAAS or KAAS-like format are supported up to February 2022).  
**Annotations format don't need to be changed from their original output** (truncated example files can be found in `KEMET/toy/` folder)  

- if genome-scale model operations are needed, optional pre-existing genome-scale models generated with CarveMe (".xml" files) that need gap-fill could be copied in the `KEMET/models/` folder, created after setup process.  

**IMPORTANT NOTE:  
KEGG annotation and GSMM files names MUST be the same of the genome they refers to, except for the extension  
(e.g. bin1.fasta, bin1.emapper.annotations, and bin1.xml)**

`KEMET/KEGG_MODULES/` folder presence as in GitHub is **necessary** for script usage. It contains files that represent KEGG Modules block structure (REF: [KEGG MODULE resource](https://www.genome.jp/kegg/module.html)); missing KO orthologs are deduced from these structures.  

Other "custom" Modules could be added to that folder, if formatted in the proper way.

3) Fill in the textual intruction file "genomes.instruction":  
excluding the header, each line should have a **tab-separated** indication of:

|MAG/Genome FASTA|Taxonomic indication|Metabolic model universe|
|---|---|---|  

- The MAG/Genome FASTA indicate the MAG/Genome of interest file name (e.g. bin1.fasta)  

- The taxonomic indication is taken from the KEGG Brite taxonomic indication (specifically from the C-level, that coincide with NCBI's phylum level most of the times) (REF: [BRITE Organism table](https://www.kegg.jp/brite/br08601)) (e.g. Actinobacteria)  

- Metabolic model universe comprehend grampos, gramneg, archaea or other custom universe (this is an optional indication needed for GSMM de-novo reconstruction)  

4) Fill in other instruction text files.  
If HMM-analyses are desired, these need either the "module_file.instruction" or the "ko_file.instruction" files as follows, depending on the desired MODE OF USE (which needs to be specified with the `--hmm_mode MODE` parameter)  

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

`--skip_hmm`: use this parameter to stop after KEGG MODULES Completeness Evaluation. The only output would be the organized tables of metabolic potential from previous annotation softwares.  

`--hmm_mode [MODE]`: if the HMM analysis is desired, this parameter is required to indicate which subset of KOs to search with profile HMMs (`onebm, module, kos`).  

`--skip_gsmm`: use this parameter to stop after HMM analysis.  

`--gsmm_mode [MODE]`: if GSMM operations with HMM-hits are desired, this parameter is required to indicate whether to perform de-novo GSMM reconstruction or gapfilling, with `[MODE]` equal to `denovo` or `existing`, respectively .  

-----
Other parameters are not required, and are described in the following paragraph:  

`-v`: print more informations, for progress reporting purposes & more info.  
`-q`: print less informations and silence MAFTT and HMMER soft-errors (suggested).  

`--as_kegg`: change the way incomplete KEGG MODULES are considered in the metabolic potential recap tables, following KEGG-Mapper convention, i.e. Modules with less than 3 blocks that have at least 1 incomplete block are marked as INCOMPLETE regardless.  

`--update_taxonomy_codes`: when downloading a new BRITE taxonomy with `setup.py`, this parameter is necessary to have an updated species (BRITE E-class) list for each C-class, that is to have an updated scope for taxonomic-informed HMM creation.  

`--threshold_value [VALUE]`: modify this to use a different quality filter to differentiate between HMM hits (default: 0.43).  

`--skip_nt_download`: avoid the KEGG GENES sequences download, e.g. if they were previously downloaded.  
`--skip_msa_and_hmmbuild`: avoid MAFFT and HMMER hmmbuild commands launch, e.g. if multi-alignment and pHMMs were done generated.  
`--retry_nhmmer`: use this parameter after a full pipeline run, to re-run nHMMER commands properly, e.g. in order to use a different scoring threshold.  

-----
# Output descriptions

## 1) KEGG Modules Completeness evaluation
Starting from previously obtained functional annotation, the MAG/Genome in input is evaluated in terms of KEGG Orthologs (KO) presence. This evaluation is performed in the framework of KEGG Modules, i.e. manually defined functional units composed of KOs, in order to recapitulate (meta)genomic potential.  
KEMET output regarding this task is composed of 2 files for each MAG/Genome of interest.  
They are located in the `KEMET/report_tsv/` and `KEMET/report_txt/` folders respectively:

### - reportKMC_FASTA.tsv  
This output is a tabular file with infomations regarding each KEGG Module, indicating the metabolic potential of the MAG/genome defined with the `FASTA` name.  

Each line includes the tab-separated info as in the following example table:  

|Module_id|Module_name|Completeness|complete/total blocks|missing KOs|KOs present|  
| --- | --- | --- | --- | --- | --- |
|M00029|Urea cycle|1 BLOCK MISSING|4__5|K01948,K14681|K00611,K01940,K01755,K01476|  

- The **Completeness** indications are accordingly: "INCOMPLETE", "2 BLOCKS MISSING", "1 BLOCK MISSING", or "COMPLETE". 

- **complete/total blocks** is indicated with the format "COMPLETE__TOTAL" (with two underscores).

### - reportKMC_FASTA.txt  
This output is a flat file with indication of KEGG MODULES completeness for every Module, up to the block level. It gives info on which sequential step of the Module path has missing KOs.

-----
## 2) HMM-based ortholog search
KEMET performs bulk nucleotidic sequences download from KEGG GENES using [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html). For license terms see [this site](https://www.kegg.jp/kegg/legal.html). API service is available free of charge to academic users only. If users prefer different download options they are encouraged to request a KEGG FTP subscription.  
Downloaded GENES sequences are filtered (all unique sequences are considered once for HMM building), aligned using MAFFT multi-sequence aligner and a profile is created using [HMMer suite](http://hmmer.org/).
The nucleotidic profiles obtained are further searched in the MAG/Genome of interest.  
As a default, a threshold value is imposed in order to enrich for complete profiles while not including hits resulting from partial sequences.  
Only hits with a score that surpass the threshold are considered proper hits, resulting in the presence of KO(s) of interest in the MAG/Genome sequences.  

Information regarding HMM hits is included in the output files:

### - FASTA_HMM_hits.txt  

A tabular file including HMM hits of a **single MAG/Genome**, defined with the `FASTA` name. It contains informations on the hits in the form:

|KO|corr_score, e-value|contig_name|strand|genome_left_bound|genome_right_bound|profile_lenght|begin_of_HMMsequence_hit|end_of_HMMsequence_hit|
|---|---|---|---|---|---|---|---|---|

- **corr_score** is a metric that describes HMM profile scoring, corrected on the sequence lenght of that profile.  

### - file_recap_DATE.tsv  

After a single KEMET run, a tabular summary file is generated. It **includes every "_HMM_hits" file information** and incorporates them in a single table.

Moreover, the file includes further fields:  

|frame|seq|xseq|
|---|---|---|  

- **frame** indicates the most likely translated reading-frame.  

- **seq** is the nucleotidic sequence as retrieved from the MAG/Genome.

- **xseq** is the translated aminoacidic sequence derived from HMM seq using the Bacterial/Archaeal translation table (t11).  

-----
## 3) Genome-scale metabolic model gapfilling
The script connects missing KOs content, retrieved via HMM hits, to reactions in the BiGG namespace (ModelSEED namespace will be added in a next release).  

Based on the `--gsmm_mode` parameter it operates in two different ways:  

`--gsmm_mode denovo` allows an automatic gene-calling from MAG/Genome sequences using Prodigal, and automatically adds the hits retrieved with HMMs to proteins multiFASTA (.faa) files.  

After that, KEMET performs a [CarveMe](https://github.com/cdanielmachado/carveme) reconstruction including these newly found sequences.  
**NOTE** The usage thus described is subject to CarveMe dependences, including the IBM CPLEX Optimizer. More regarding the dependencies can be read about CarveMe installation procedure [here](https://carveme.readthedocs.io/en/latest/installation.html).

Using this mode, the newly generated gene prediction and GSMM are included in the `KEMET/de_novo_models` folder.  

`--gsmm_mode existing` allows the identified reactions to be incorporated in esisting genome-scale metabolic models (GSMMs) previously generated with CarveMe, if those are missing.  

At the moment (February 2022) the only tested way to add reaction to pre-existing GSMMs is via the [ReFramed](https://github.com/cdanielmachado/reframed) package.
Further improvement would permit adding it through the [cobrapy](https://github.com/opencobra/cobrapy) platform.  

Informations regarding reaction gapfilling (if performed using the `--gsmm_mode existing` parameter) are included in several output files:

### - bigg_log_FASTA.txt
A flat-file with the indication of every BiGG reaction that potentially could be added to the model in input, defined with the `FASTA` name. The BiGG reactions are included in a one per line format.  

### - FASTA_added_reactions
A flat-file with the reactions that were actually added for a given MAG/Genome-derived GSMM, defined with the `FASTA` name. Reaction names are indicated one per line. followed by the respective reaction string.  

### - gapfilled model
Individual GSMMs are saved again after the gapfilling procedure with new reactions and metabolites content as `FASTA_KEGGadd_DATE.xml`, where `FASTA` follows the input definition and `DATE` includes the day of analysis. Files generated this way are stored in the `KEMET/model_gapfilled/` folder.  

-----
# Credits

Developed by Matteo Palù at Università degli Studi di Padova (2020-2022).

## _help pages_

### _setup.py help_
```
usage: setup.py [-h] [-k] [-u] [-G]

Setup command for KEMET pipeline. Create folders and generate/update KEGG
Module .kk database

optional arguments:
  -h, --help           show this help message and exit
  -k, --set_kk_DB      Choose this option in order to create KEGG Module DB
                       (.kk files), in order to perform KEGG Modules
                       COmpleteness evaluation.
  -u, --update_kk_DB   Choose this option in order to update already existing
                       KEGG Module DB (.kk files).
  -G, --gapfill_usage  Choose this option in order to create and download
                       required folders for Gapfilling, follow-up of HMM
                       search.
```

### _kemet.py help_
```
usage: kemet.py [-h] -a {eggnog,kaas,kofamkoala} [--update_taxonomy_codes]
                [-I PATH_INPUT] [-k] [--skip_hmm]
                [--hmm_mode {onebm,modules,kos}]
                [--threshold_value THRESHOLD_VALUE] [--skip_nt_download]
                [--skip_msa_and_hmmbuild] [--retry_nhmmer] [--skip_gsmm]
                [--gsmm_mode {existing,denovo}] [-O PATH_OUTPUT] [-v] [-q]
                [--log]
                FASTA_file

    KEMET pipeline:
    1) Evaluate KEGG Modules Completeness for given genomes.
    2) HMM-based check for ortholog genes (KO) of interest after KEGG Module Completeness evaluation.
    3) Genome-scale model gapfill with nucleotidic HMM-derived evidence, for KOs of interest.
    

positional arguments:
  FASTA_file            Genome/MAG FASTA file as indicated in the "genomes.instruction" -
                        points to files (in "KEGG_annotations") comprising KO annotations, associated with each gene.

optional arguments:
  -h, --help            show this help message and exit
  -a {eggnog,kaas,kofamkoala}, --annotation_format {eggnog,kaas,kofamkoala}
                        Format of KO_list.
                        eggnog: 1 gene | many possible annotations;
                        kaas: 1 gene | 1 annotation at most;
                        kofamkoala: 1 gene | many possible annotations
  --update_taxonomy_codes
                        Update taxonomy filter codes - WHEN TO USE: after downloading a new BRITE taxonomy with "setup.py".
  -I PATH_INPUT, --path_input PATH_INPUT
                        Absolute path to input file(s) FOLDER.
  -k, --as_kegg         Return KEGG-Mapper output for the Module Completeness evaluation.
  --skip_hmm            Skip HMM-driven search for KOs & stop after KEGG Modules Completeness evaluation.
  --hmm_mode {onebm,modules,kos}
                        Choose the subset of KOs of interest for HMM-based check.
                        By default, the KOs already present in the functional annotation are not checked further.
                        
                        onebm: search for KOs from KEGG Modules missing 1 block;
                        modules: search for KOs from the KEGG Modules indicated in the "module_file.instruction" file, 1 per line
                            (e.g. Mxxxxx);
                        kos: search for KOs indicated in the "ko_file.instruction" file, 1 per line
                            (e.g. Kxxxxx)
                                                
  --threshold_value THRESHOLD_VALUE
                        Define a threshold for the corrected score resulting from HMM-hits, which is indicative of good quality.
  --skip_nt_download    Skip downloading KEGG KOs nt sequences.
  --skip_msa_and_hmmbuild
                        Skip MAFFT and HMMER hmmbuild commands.
  --retry_nhmmer        Move HMM-files and re-run nHMMER command.
  --skip_gsmm           Skip GSMM operations, gapfill or de-novo model creation, & stop after HMM-driven search for KOs.
  --gsmm_mode {existing,denovo}
                        Choose the methods of GSMM operation.
                            (This method won't be performed if "--hmm_mode kos" was chosen)
                        existing: use pre-existing CarveMe GSMM to add reactions content connected to HMM-derived KOs;
                        denovo: generate a new CarveMe GSMM, performing gene prediction and adding HMM-derived hits from the chosen HMM-mode.
                                                
  -O PATH_OUTPUT, --path_output PATH_OUTPUT
                        Absolute path to ouput file(s) FOLDER.
  -v, --verbose         Print more informations - for debug and progress.
  -q, --quiet           Silence soft-errors (for MAFFT and HMMER commands).
  --log                 Store KEMET commands and progress during the execution in a log file.
```
