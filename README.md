# KEMET
KEgg Module Evaluation Tool

## General Use

0) install dependences

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
```

1) Clone the git repository.

```
git clone https://github.com/Matteopaluh/KEMET.git
```

2) Run the _setup.py_ script once, specifying the scope (if HMM analysis is needed, add "-H" option; for genome-scale models possibilities, add "-G" option).
This will set the folders in which different input and outputs will be stored.

The "KEGG_MODULES" directory from the Git is mandatory for script usage.

3) Setup input into proper paths:
- KEGG KOs annotations (derived from different sources) should be copied in the "KEGG_mappings" folder, which is created in the setup process.
- The script needs an indication of KEGG mappings format (eggNOG, KofamKOALA, and KAAS or KAAS-like supported up to now 26/05/21).
- Two tipes of output are possible, alone or in combination, as specified later.

4) The scripts (_kmc.py_, _hmm.py_, _gsmm.py_) are meant to be utilized in a sequential order: latter step scripts need previous steps outputs.

### Single Scripts Use

### 1) kmc.py

"KEMET/KEGG_MODULES/kk_files" contains ".kk" files that represent the block structure of Modules; those files are scanned by the script in order to identify the missing KO orthologs of those structures.

- KEGG KOs annotations (derived from different sources) should be copied in the "KEMET/KEGG_mappings/" folder, which is created after setup process.

- Other "custom" Modules could be added to that folder, with a proper format. (REF: [KEGG MODULE resource](https://www.genome.jp/kegg/module.html))

**Batch use**:
for f in $(ls ./KEGG_mappings/); do ./kmc.py $f <OPTIONS>; done

#### Outputs
##### .txt (option: -o txt)
A flat file with indication of KMC for every Module, up to the block level. It gives info on which sequential step of the Module path has missing KOs.
##### .tsv (option: -o tsv)
A tab-separated file. Infos for each line are the following:

Module_id; Module_name; Status (i.e. (IN-)COMPLETE/1-2 BLOCK MISSING); complete/total blocks (C__T); KOs missing; KOs present

##### .txt and .tsv together (option: -o txt+tsv)
Generate both of previous outputs.


### 2) hmm.py

The script performs bulk nt sequences download using KEGG API.
Afterwards those sequences are filtered, aligned using a multi-sequence aligner and a profile is created using [HMMer suite](http://hmmer.org/).
The nucleotidic profiles obtained are further searched in the MAG/Genome of interest.
Only hits satisying given criteria are considered proper hits, signifying the presence of a given KO in the MAG/Genome genomic sequence, and included in the outputs.

- Compile "genomes.instruction" tabular file (created with "setup.py") with the ID of MAG/Genome of interest, the KEGG Brite taxonomic indication (C-level, that coincide with NCBI's phylum level most of the times), and eventually an indication of the metabolic-model universe (grampos, gramneg, archaea or such). **This can be done using helping script "add-genomes-info.py"**

- Decide MODE OF USE, checking between different group of orthologs:
    1. KOs from KEGG Modules missing 1 block								(OPTION: --onebm_modules_list)
    2. KOs from a fixed list of KEGG Modules, indicated one per line in "module_file.instruction"	(OPTION: --fixed_modules_list)
    3. KOs from a fixed list of orthologs, indicated one per line in "ko_file.instruction"		(OPTION: --fixed_ko_list)

- Eventually, a different threshold value can be set (OPTION: --threshold_value)

- In case of pre-downloaded nucleotidic sequences // pre-computed alignment and profile creation, different steps of the program can be skipped (OPTIONS: --skip_nt_download, --skip_msa_and_hmmbuild)

- It is possible to re-analyse data in the final part of the program, using different scoring parameters (OPTION: --retry_nhmmer)

#### Outputs
##### MAG_HMM_hits.txt files
A tab-separated file originated from hits of a single MAG/Genome. It gives informations on the hits regarding:

KO; score (corrected by profile lenght), e-value; contig_name; strand; genome_left_bound; genome_right_bound; profile_lenght; begin_of_HMMsequence_hit; end_of_HMMsequence_hit

##### file_recap_DATE.tsv
General tab-separated summary file. Includes every "_HMM_hits" file information and group them;
moreover it includes a field for the most likely translated reading-frame, the nucleotidic sequence as retrieved from the MAG/Genome, as well as the translated aminoacidic sequence using Bacterial/Archaeal translation.

**Batch use**:
for f in $(ls ./genomes/); do f1="${f##*/}"; f2=${f1%.*}; ./hmm.py $f2 OPTIONS; done


### 3) gsmm.py

The script connects missing KO content, identified via HMM-hits, to reactions in the BiGG namespace.
Furthermore it adds those reactions to Genome-scale metabolic models (GSMMs) generated with CarveMe, if missing.

At the moment (31/03/21) the only tested way to add reaction is via the [ReFramed](https://github.com/cdanielmachado/reframed) package.
Further improvement would permit adding it through the [cobrapy](https://github.com/opencobra/cobrapy) platform.

- CarveMe GSMMs (".xml" files) should be copied in the "KMC/models/" folder, which is created after setup process.

- The same "genomes.instruction" tabular file as "kmc-hmm.py" is utilized, to get the ID of MAG/Genome of interest.

- Decide MODE OF USE, checking between different group of orthologs:
    1. KOs from KEGG Modules missing 1 block								(OPTION: --onebm_modules_list)
    2. KOs from a fixed list of KEGG Modules, indicated one per line in "module_file.instruction"	(OPTION: --fixed_modules_list)

#### Outputs
##### bigg_log_MAG.txt files
A flat-file with the indication of every BiGG reaction that could be added to the model in input.
The BiGG reactions are indicated one per line.

##### MAG_added_reactions.txt files
A flat-file with the actual added reactions for a given MAG/Genome-derived GSMM.
The name of the reaction is followed by the reaction string.

##### gapfilled models files
Individual GSMMs are saved again with new reaction and metabolite content as "MAG_KEGGadd_DATE.xml".


## Credits

Developed by Matteo Palù at Università degli Studi di Padova (2020-2021).

### _help pages_

#### _setup.py help_
```
usage: setup.py [-h] [-k] [-u] [-H] [-G] [-v]

Setup command for KMC+HMM analysis. Create folders and generate/update KEGG
Module .kk database

optional arguments:
  -h, --help           show this help message and exit
  -k, --set_kk_DB      Choose this option in order to create KEGG Module DB
                       (.kk files), in order to use KMC.
  -u, --update_kk_DB   Choose this option in order to update already existing
                       KEGG Module DB (.kk files).
  -H, --hmm_usage      Choose this option in order to create and download
                       required folders for HMM analysis, follow-up of KMC.
  -G, --gapfill_usage  Choose this option in order to create and download
                       required folders for Gapfilling, follow-up of KMC+HMM.
  -v, --verbose        Print more informations - for debug or log.
```

#### _kmc.py help_
```
usage: kmc.py [-h] -a {eggnog,kaas,kofamkoala} [-I PATH_INPUT]
              [-o {txt,tsv,txt+tsv}] [-k] [-O PATH_OUTPUT] [-v]
              KOfile

Evaluate KEGG Modules Completeness for given genomes.

positional arguments:
  KOfile                File comprising KO annotations, associated with each
                        gene.

optional arguments:
  -h, --help            show this help message and exit
  -a {eggnog,kaas,kofamkoala}, --annotation_format {eggnog,kaas,kofamkoala}
                        Format of KO_list. eggnog: 1 gene | many possible
                        annotations; kaas: 1 gene | 1 annotation at most.
  -I PATH_INPUT, --path_input PATH_INPUT
                        Absolute path to input file(s) FOLDER.
  -o {txt,tsv,txt+tsv}, --output {txt,tsv,txt+tsv}
                        Output format for KMC summary. txt: more level-
                        detailed, worse in recap; tsv: best at recap, easily
                        parsable for downstream analysis; txt+tsv: both of the
                        above.
  -k, --as_kegg         Return KEGG Mapper output for the Module completeness.
  -O PATH_OUTPUT, --path_output PATH_OUTPUT
                        Absolute path to ouput file(s) FOLDER.
  -v, --verbose         Print more informations - for debug and progress.
```

#### _hmm.py help_
```
usage: hmm.py [-h] [--update_taxonomy_codes] [--onebm_modules_list]
              [--fixed_modules_list] [--fixed_ko_list]
              [--threshold_value THRESHOLD_VALUE] [--skip_nt_download]
              [--do_aa_download] [--skip_msa_and_hmmbuild] [--retry_nhmmer]
              [-I PATH_INPUT] [-v] [-q]
              MAG_genome_FASTA

HMM-based check for ortholog genes after KEGG Module Completeness evaluation.

positional arguments:
  MAG_genome_FASTA      Run HMM-based search for KOs in Modules of interest in
                        the genome indicated with this expression.

optional arguments:
  -h, --help            show this help message and exit
  --update_taxonomy_codes
                        Update taxonomy filter codes - WHEN TO USE: after
                        downloading a new BRITE taxonomy with "setup.py".
  --onebm_modules_list  Use all KEGG Modules missing 1 block for the HMM-based
                        check (only KOs missing in the indicated MAG/Genome).
  --fixed_modules_list  Use a fixed list of KEGG Modules to use for the HMM-
                        based check (only KOs missing in the indicated
                        MAG/Genome).
  --fixed_ko_list       Use a fixed list of KOs to use for the HMM-based check
                        (only KOs missing in the indicated MAG/Genome).
  --threshold_value THRESHOLD_VALUE
                        Define a threshold for the corrected score resulting
                        from HMM-hits, which is indicative of good quality.
  --skip_nt_download    Skip downloading KEGG KOs nt sequences.
  --do_aa_download      Allow downloading KEGG KOs aa sequences.
  --skip_msa_and_hmmbuild
                        Skip MAFFT and HMMER hmmbuild commands.
  --retry_nhmmer        Move HMM-files and re-run nHMMER command.
  -I PATH_INPUT, --path_input PATH_INPUT
                        Absolute path to input file(s) FOLDER.
  -v, --verbose         Print more informations - useful for debug and
                        progress.
  -q, --quiet           Silence soft-errors (for MAFFT and HMMER commands).
```

#### _gsmm.py help_
```
usage: gsmm.py [-h] -r {reframed,cobrapy} -d {bigg,modelseed} [--de_novo_GSMM]
               [--onebm_modules_list] [--fixed_modules_list] [-I PATH_INPUT]
               [-O PATH_OUTPUT] [--log] [-v]
               MAG_genome_FASTA

Genome-scale model reaction-addition after KEGG Modules Completeness
evaluation and HMM-evidence for given genomes.

positional arguments:
  MAG_genome_FASTA      Use the fasta ID as indicated .

optional arguments:
  -h, --help            show this help message and exit
  -r {reframed,cobrapy}, --reaction_addition_method {reframed,cobrapy}
                        Selection of model I/O program for reaction addition.
  -d {bigg,modelseed}, --ontology_database_selection {bigg,modelseed}
                        Selection of reaction ontology for gap-fill addition
                        in GSMM.
  --de_novo_GSMM        Use all KEGG Modules missing 1 block for the HMM-based
                        check (only KOs missing in the indicated MAG/Genome).
  --onebm_modules_list  Use all KEGG Modules missing 1 block for the HMM-based
                        check (only KOs missing in the indicated MAG/Genome).
  --fixed_modules_list  Use a fixed list of KEGG Modules to use for the HMM-
                        based check (only KOs missing in the indicated
                        MAG/Genome).
  -I PATH_INPUT, --path_input PATH_INPUT
                        Absolute path to input file(s) FOLDER.
  -O PATH_OUTPUT, --path_output PATH_OUTPUT
                        Absolute path to ouput file(s) FOLDER.
  --log                 Store reaction-addition commands / info as plain text
                        files.
  -v, --verbose         Print more informations - for debug or log.
```
