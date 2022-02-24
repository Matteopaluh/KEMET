# Examples of KEGG annotations

The files included in this folder are examples of input KEGG annotations, and this README includes a discussion about the different formats, along with info on how to use the different output files within KEMET. Said files only include the first few lines of what an actual functional annotation file would be made of, only to show their formatting differences.  

The files extension **does not need to be modified**, neither the rest of file name, except it does not refer to the same genome in input (e.g. bin1.fasta/.fa/.fna needs to be selected).  

-----
- `bin1.emapper.annotations` is an example of one of the files that can be generated using either eggNOG mapper command line program or the [web server](http://eggnog-mapper.embl.de/) (the other being, in this case, `bin1.seed_orthologs`, not used by KEMET).  
To utilize this type of annotation within `kemet.py` script, the user will need to add the `-a` or `--annotation_format` argument, then write "eggnog": `-a eggnog`  

- `bin1.ko` is an example of a file generated via KAAS (KEGG Automatic Annotation Server), which can be accessed [here](https://www.genome.jp/tools/kaas/).  
To utilize this type of annotation within `kemet.py` script, the user will need to add the `-a` or `--annotation_format` argument, then write "kaas": `-a kaas`  

- `bin1.txt` is an example of a file generated via KofamKOALA, which is accessed either via command line or via [web server](https://www.genome.jp/tools/kofamkoala/).  
Specifically, if the request is executed via web server, this file is the output that can be downloaded with the **"result file"** option.  
To utilize this type of annotation within `kemet.py` script, the user will need to add the `-a` or `--annotation_format` argument, then write "kofamkoala": `-a kofamkoala`  

- `bin1_ko.txt` is an example of a KofamKOALA output downloaded with the **"Input file for KEGG Mapper"** option.  
This file has the same format of those generated via KAAS, i.e. a tabular file with two fields, one indicating a gene-id, the other either being empty or indicating either a KEGG Ortholog (KO). Within the KEMET vocabulary, this is referred to as "kaas-like".  
To utilize this type of annotation within `kemet.py` script, the user will need to add the `-a` or `--annotation_format` argument, then write "kaas": `-a kaas`