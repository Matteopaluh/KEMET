# KEMET
KEgg Module Evaluation Tool

## Script usage description

The _kemet.py_ script works in a sequential order for each genome utilized. It is well suited for single genomes usage as well as metagenome-spanning analyses.
KEMET works as a command line tool that serves three main functions:

1) Evaluate KEGG Modules Completeness, producing organized tables for each MAG/Genome of interest. 
2) HMM-based searches for ortholog genes (KO) of interest chosen from a subset after KEGG Module Completeness evaluation.
3) Genome-scale model gapfill with evidence from nucleotidic HMM searches, regarding KOs of interest.

## General Use - Conda environment

0) Install dependences

In a conda environment with python >= 3.6, use the following installation commands:

```
conda install -c bioconda mafft
sudo apt install hmmer
pip install reframed
pip install cherrypy
```

OPTIONAL
```
pip install carveme=1.4.1
conda install -c bioconda diamond
conda install -c bioconda prodigal
```

1) Clone the git repository.

```
git clone https://github.com/Matteopaluh/KEMET.git
```

2) Enable execution of the scripts.
```
chmod +x .
```

3) Run the _setup.py_ script, specifying the scope (for genome-scale models possibilities, add the "-G" option).
This will set the folders in which different input and outputs will be stored.

4) Setup input files into proper paths:
- MAG/Genome sequences should be copied in the "/genomes/" folder, which is created after the setup process.
- KEGG KOs annotations (derived from different sources) should be copied in the "/KEGG_annotations/" folder, which is created in the setup process.
- if genome-scale model operations are needed, CarveMe genome-scale models (".xml" files) to gap-fill could be copied in the "KEMET/models/" folder, which is created after setup process.
**IMPORTANT NOTE: each KEGG annotation and GSMM files should be called as the genome it refers to, except for the extension (e.g. bin1.fasta, bin1.emapper.annotations, and bin1.xml)**
- The script needs an indication of KEGG annotation format (eggNOG, KofamKOALA, and KAAS or KAAS-like supported up to now 07/10/21).

The "KEMET/KEGG_MODULES/kk_files/" directory presence from the Git is **mandatory** for script usage. It contains ".kk" files that represent the block structure of Modules (REF: [KEGG MODULE resource](https://www.genome.jp/kegg/module.html)); these files are "scanned" in order to identify the missing KO orthologs of those structures.

- Other "custom" Modules could be added to that folder, with a proper format.

5) Fill in the intruction text file "genomes.instruction"; excluding the header, each line should have a tab-separated indication of:
- MAG/Genome FASTA file name (e.g. bin1.fasta)
- the KEGG Brite taxonomic indication (C-level, that coincide with NCBI's phylum level most of the times) (REF: [BRITE Organism table](https://www.kegg.jp/brite/br08601)) (e.g. Actinobacteria)
- an indication of metabolic-model universe (grampos, gramneg, archaea or such)(optional - needed for GSMM de-novo reconstruction)

6) Fill in the other instruction text files. If HMM-analyses are desired, either the "module_file.instruction" or the "ko_file.instruction" files, depending on the desired MODE OF USE (which need to be specified in the --hmm_mode parameter):

    1. KOs from KEGG Modules missing 1 block								(No need to fill instruction files)
    2. KOs from a fixed list of KEGG Modules	              (Need to be indicated one per line in the "module_file.instruction" file)
    3. KOs from a fixed list of orthologs                   (Need to be indicated one per line in the "ko_file.instruction" file)

7) Launch the _kemet.py_ script as a command-line with the arguments of choice.

## Command line (minimal, only required arguments)
```
./kemet.py FASTA_file -a ANNOTATION_FORMAT --skip_hmm --hmm_mode MODE --skip_gsmm --gsmm_mode MODE
```

**FASTA_file**: indication of the MAG/Genome of interest FASTA (with or without path indication e.g. genomes/bin1.fasta)\
**-a ANNOTATION_FORMAT**: indication of the KEGG annotation format (either eggnog/kaas/kofamkoala), which is mandatory to generate KEGG MODULES recap tables.\
**--skip_hmm**: use this parameter to stop after KEGG MODULES Completeness Evaluation, to have only the tables from previous annotation softwares.\
**--hmm_mode MODE**: if the HMM analysis is desired, this parameter is required to indicate which subset of KOs to search with profile HMMs.\
**--skip_gsmm**: use this parameter to stop after HMM analysis.\
**--gsmm_mode MODE**: if GSMM operation with HMM-hits are desired, this parameter is required to indicate whether to perform de-novo GSMM reconstruction or gapfilling.\

Other parameters are not required, but are described in the following paragraph:\
**VERBOSITY**\
**-v**: print more informations, for progress reporting purposes & more info.\
**-q**: print less informations and silence MAFTT and HMMER soft-errors.\

**--as_kegg**: change the way incomplete KEGG MODULES are considered in the recap tables, following KEGG-Mapper convention, i.e. Modules with less than 3 blocks that have at least 1 incomplete block are marked as INCOMPLETE regardless.\
**--update_taxonomy_codes**: when downloading a new BRITE taxonomy with _setup.py_, this parameter is necessary to have an updated species (BRITE E-class) list for each C-class, that is to have an updated scope for taxonomic-informed HMM creation.\
**--threshold_value VALUE**: use this parameter to have a different HMM corrected score threshold to differentiate between proper hits (default: 0.4).\
**--skip_nt_download**: skip downloading KEGG KOs nt sequences, e.g. if they were previously downloaded.\
**--skip_msa_and_hmmbuild**: skip MAFFT and HMMER hmmbuild commands, e.g. if multi-alignment and profile HMMs were done previously.\
**--retry_nhmmer**: use this parameter after a full pipeline run, to re-run nHMMER commands properly, e.g. in order to use a different threshold.\

## Output descriptions

### 1) KEGG Modules Completeness evaluation
An automated search is performed for each KEGG Module, i.e. manually defined functional units composed of KEGG orthologs (KOs), starting from previously performed functional annotation.
KEMET outputs 2 types of files for each MAG/Genome of interest, located in the "/report_tsv/" and "/report_txt/" folders respectively:

- **"reportKMC_"+FASTA+".tsv"**: a tab-separated file with infomations about each KEGG Module. Each line has the following info:
Module_id; Module_name; Completeness (i.e. (IN-)COMPLETE/1-2 BLOCK MISSING); complete/total blocks (COMPLETE__TOTAL); KOs missing; KOs present

- **"reportKMC_"+FASTA+".txt"**: a flat file with indication of KEGG MODULES completeness for every Module, up to the block level. It gives info on which sequential step of the Module path has missing KOs.

### 2) HMM-based ortholog search
The script performs bulk nucleotidic sequences download using [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html).
Those sequences are filtered (all unique sequences are considered once for HMM building), aligned using a multi-sequence aligner and a profile is created using [HMMer suite](http://hmmer.org/).
The nucleotidic profiles obtained are further searched in the MAG/Genome of interest.
As a default, a threshold value is imposed in order to enrich for more complete profiles instead of partial sequences.
Only hits with a score over threshold are considered proper hits, resulting in the presence of KO(s) of interest in the MAG/Genome sequences.

This information is included in the output:

- **FASTA+"_HMM_hits.txt"**: a tab-separated file originated from HMM hits of a single MAG/Genome. It contains informations on the hits regarding:
KO; score (corrected by profile lenght), e-value; contig_name; strand; genome_left_bound; genome_right_bound; profile_lenght; begin_of_HMMsequence_hit; end_of_HMMsequence_hit

- **"file_recap_"+DATE.tsv**: After a single KEMET run, a tab-separated summary file is generated. It includes every "_HMM_hits" file information and group them;
moreover it includes other fields:
the most likely translated reading-frame; the nucleotidic sequence as retrieved from the MAG/Genome; the translated aminoacidic sequence using Bacterial/Archaeal translation.

### 3) Genome-scale metabolic model gapfilling
The script connects missing KOs content, identified via HMM hits, to reactions in the BiGG namespace (ModelSEED namespace will be added in a next release).
Furthermore it adds those reactions to genome-scale metabolic models (GSMMs) generated with CarveMe, if missing.
Otherwise it can perform MAG/Genome gene-calling and automatically add translated sequences to multiFASTA (.faa) files. After that, the script perform a [CarveMe](https://github.com/cdanielmachado/carveme) reconstruction including these newly found sequences.\
**NOTE** This usage needs access to CarveMe dependences, including IBM CPLEX Optimizer. More regarding the dependencies can be read [here](https://carveme.readthedocs.io/en/latest/installation.html).

At the moment (07/10/21) the only tested way to add reaction to pre-existing GSMMs is via the [ReFramed](https://github.com/cdanielmachado/reframed) package.
Further improvement would permit adding it through the [cobrapy](https://github.com/opencobra/cobrapy) platform.
Informations regarding reaction gapfilling (if performed using the "--gsmm_mode existing" parameter) are included in several output files:

**bigg_log_+FASTA.txt**: a flat-file with the indication of every BiGG reaction that could be added to the model in input. The BiGG reactions are indicated one per line.

**FASTA+_added_reactions.txt**: a flat-file with the reactions actually added for a given MAG/Genome-derived GSMM. The name of the reaction is followed by the reaction string.

**gapfilled models**: individual GSMMs are saved again with new reaction and metabolite content as "FASTA+_KEGGadd_+DATE.xml". These files are located in the "/model_gapfilled/" folder.

When using the "--gsmm_mode denovo" parameter, the newly generated gene prediction and GSMM are included in the "KEMET/de_novo_models" folder.


# Credits

Developed by Matteo Palù at Università degli Studi di Padova (2020-2021).

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
usage: kemet.py [-h] -a {eggnog,kaas,kofamkoala}
                     [--update_taxonomy_codes] [-I PATH_INPUT] [-k]
                     [--skip_hmm] [--hmm_mode {onebm,modules,kos}]
                     [--threshold_value THRESHOLD_VALUE] [--skip_nt_download]
                     [--skip_msa_and_hmmbuild] [--retry_nhmmer] [--skip_gsmm]
                     [--gsmm_mode {existing,denovo}] [-O PATH_OUTPUT] [-v]
                     [-q]
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
```
