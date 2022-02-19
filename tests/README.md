# Quality tests for HMM results

## Dataset definition

### MAGs
The original dataset included several MAGs from Bioproject [PRJNA602310](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602310), mentioned in [this study](https://pubmed.ncbi.nlm.nih.gov/33714811/).

The MAGs thus indicated share the fact that their exact taxonomy is unknown.
Indeed, among the 12 selected MAGs:

- 1 is known only at *CLASS* level
- 2 are known only at *ORDER* level
- 6 are known only at *FAMILY* level
- 1 is known only at *GENUS* level
- 2 are known only at *SPECIES* level (one of them being a _Candidatus_ species)


### Complete Genomes
Moreover the dataset comprised 5 Complete Genomes:

- Escherichia coli str. K-12 substr. MG1655, complete genome [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)
- Mycoplasma genitalium G37, complete sequence [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_000908.2)
- Bacillus subtilis subsp. subtilis str. 168 complete genome [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_000964.3)
- Hungateiclostridium thermocellum ATCC 27405, complete sequence [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_009012.1)
- Shewanella oneidensis MR-1 chromosome [link](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP053946.1)

**A change was made so that the FASTA header of those genomes was changed, keeping only the indication of the Nucleotide entry**

-----
## Tests

#### **Base results**

The tests were run in a conda environment with the following versions of required packages:

```
python	    3.6.13
carveme	    1.4.1
diamond	    2.0.9
mafft	    7.475
prodigal	2.6.3
reframed	1.2.1
hmmer	    3.1b2
```

The following Prodigal command was used:

```
prodigal -i $f.fna -a ./prodigal/$f.faa -o ./prodigal/$f.gff -d ./prodigal/$f.fna -q -m -f gff
```

All initial Prodigal predictions are included in the "base_results/prodigal/" folders, that can be access after the extraction from the `test_MAG_results.tar.gz` and `test_complete_genomes_results.tar.gz` archives.

eggNOG-mapper v2 web-version, using default values, was used for functional annotation. The files are located in "base_results/eggnog/" folders, that can be accessed after the extraction from the `test_results.tar.gz` and `test_complete_genomes_results.tar.gz` archives.

-----
### **Multiple ORF removal from MAGs**

#### **ORF removal (70% and 30%, respectively bigg and small deletion) in Complete Modules**

This part was conducted to obtain "mock" MAGs similar to the original except for the deleted part.  
It therefore resulted in contig files (.fna) lacking either 30% or 70% of the sequence, starting from the first nt of the identified ORF.

#### **Multiple ORF removal from MAG ANNOTATIONS**

Removal of the same genes indicated in the [details](https://github.com/Matteopaluh/KEMET/tree/revision/tests/test_details.md), from eggnog-mapper annotations. These genes were also chosen due to their annotation being only present in a single gene.    

In the "report_tsv" files the KOs of the selected genes was changed from present to absent.  

in the "report_txt" files the KOs of the selected genes were added among the missing genes.  

The genes were removed using functions (`rewrite_MAG_removing_genes()` among others) contained in the `confrontation_results.ipynb` notebook. This code was also used to identify the boundaries indicated in the detailed tests.

The results of HMM-hits summaries were scrutinized to compare if the genes removed from annotations were predicted as present in the MAGs (both original and mock-MAGs, lacking part of the coding sequenes).  

MAGs were grouped together according to the status of the coding sequences under test:  

- The "control" group (no removal was performed on those MAGs' contigs) is is referred to with the ID number "0".
- The "big deletion" (70% ORF lenght) group is referred to with the ID number "1".
- The "small deletion" (30% ORF lenght) group is referred to with the ID number "2".  

These files can be extracted from the `test_MAG_results.tar.gz` archive.

-----
### **Multiple ORF removal from Complete genomes**

#### **ORF removal (70% and 30%, respectively bigg and small deletion) in Complete Modules**

Likewise MAGs, the same protocol was followed for complete genomes, therefore the reader is referred to the earlier paragraph.  

#### **Multiple ORF removal from Complete Genomes ANNOTATIONS**

Removal of the same genes indicated in the [details](https://github.com/Matteopaluh/KEMET/tree/revision/tests/test_details.md), from eggnog-mapper annotations. These genes were also chosen due to their annotation being only present in a single gene.    

In the "report_tsv" files the KOs of the selected genes was changed from present to absent.  

in the "report_txt" files the KOs of the selected genes were added among the missing genes.

The genes were removed using functions (`rewrite_Genome_removing_genes()` among others) written in the `confrontation_results.ipynb` notebook. This code was also used to identify the boundaries indicated in the detailed tests.  

The results of HMM-hits summaries were scrutinized to compare if the genes removed from annotations were predicted as present in the MAGs (both original and mock-genomes, lacking part of the coding sequenes).  

Complete genomes were grouped together according to the status of the coding sequences under test:

- The "control" group (no removal was performed on those complete genomes' contigs) is is referred to with the ID number "4".
- The "big deletion" (70% ORF lenght) group is referred to with the ID number "5".
- The "small deletion" (30% ORF lenght) group is referred to with the ID number "6".  

These files can be extracted from the `test_complete_genomes_results.tar.gz` archive.