# Quality tests for KMC-HMM script

## Dataset definition

##### MAGs
The original dataset included several MAGs from Bioproject [PRJNA602310](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602310), mentioned in [this study](https://pubmed.ncbi.nlm.nih.gov/33714811/).

The MAGs thus indicated share the fact that their exact taxonomy is unknown.
Indeed, among the 12 selected MAGs:

- 1 is known only at *CLASS* level
- 2 are known only at *ORDER* level
- 6 are known only at *FAMILY* level
- 1 is known only at *GENUS* level
- 2 are known only at *SPECIES* level (one of them being a _Candidatus_ species)


##### Complete Genomes
Moreover the dataset comprised 5 Complete Genomes:

- Escherichia coli str. K-12 substr. MG1655, complete genome [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)
- Mycoplasma genitalium G37, complete sequence [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_000908.2)
- Bacillus subtilis subsp. subtilis str. 168 complete genome [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_000964.3)
- Hungateiclostridium thermocellum ATCC 27405, complete sequence [link](https://www.ncbi.nlm.nih.gov/nuccore/NC_009012.1)
- Shewanella oneidensis MR-1 chromosome [link](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP053946.1)

**A change was made so that the FASTA header of those genomes was changed, keeping only the indication of the Nucleotide entry**

## Tests

#### Base results

From MAGs data, the following were chosen:

METABAT_AS10tlH2TH_158,
METABAT_AS26fmACSIPLY_35,
METABAT_AS15tlH2ME_193,
METABAT_AS07pgkLD_51,
METABAT_AS15tlH2ME_52,
METABAT_AS07pgkLD_55,
METABAT_AS15tlH2ME_127,
METABAT_AS24abBPME_148,
METABAT_AS07pgkLD_225,
METABAT_AS05jafATM_34,
METABAT_AS20ysBPTH_14,
METABAT_AS20ysBPTH_159

The tests were run in a conda environment with the following versions of required packages:

```
python	3.6.13
carveme	1.4.1
cherrypy	18.6.0
diamond	2.0.9
mafft	7.475
prodigal	2.6.3
reframed	1.2.1
hmmer	3.1b2
```

The following Prodigal command was used:

```
prodigal -i $f.fna -a ./prodigal/$f.faa -o ./prodigal/$f.gff -d ./prodigal/$f.fna -q -m -f gff
```

All initial Prodigal predictions are included in the "base_results/prodigal/" folders.

eggNOG-mapper web-version, using default values was used for functional annotation. The files are located in "base_results/eggnog/" folders.

#### Multiple ORF removal from MAGs:

##### ORF removal (70% and 30%, respectively bigg and small deletion) in Complete Modules

This part was conducted to obtain "mock" MAGs similar to the original except for the deleted part;
it therefore resulted in contig files (.fna) lacking either 30% or 70% of the sequence, starting from the first nt of the ORF indicated in the following bullet points.

1. METABAT_AS10tlH2TH_158:

1.1. K00036 glucose-6-phosphate 1-dehydrogenase [EC:1.1.1.49 1.1.1.363]
- present in M00004,M00006,M00008
- present in gene 115385_AS10_1
    - in MAGs contigs: >115385_AS10 (MINUS_STRAND) - genomic boundaries: 2-1331
- sequence present in Prodigal data
- pre-test: M00008 COMPLETE 4/4; M00004,M00006 1 BLOCK MISSING 7/8 2/3)

1.2. K00620 glutamate N-acetyltransferase / amino-acid N-acetyltransferase [EC:2.3.1.35 2.3.1.1]
- present in M00028
- present in gene 72067_AS10_35
    - in MAGs contigs: >72067_AS10 (MINUS_STRAND) - genomic boundaries: 41277-42495
- sequence present in Prodigal data
- pre-test: M00028 COMPLETE 5/5

1.3. K03340 diaminopimelate dehydrogenase [EC:1.4.1.16]
- present in M00526
- present in gene 174723_AS10_25
    - in MAGs contigs: >174723_AS10 (MINUS_STRAND) - genomic boundaries: 24849-25836
- sequence present in Prodigal data
- pre-test: M00526 COMPLETE 6/6

2. METABAT_AS26fmACSIPLY_35:

2.1. K00850 6-phosphofructokinase 1 [EC:2.7.1.11]
- present in M00001,M00345
- present in gene 20937_AS26_26
    - in MAGs contigs: >20937_AS26 (MINUS_STRAND) - genomic boundaries: 29040-30015
- sequence present in Prodigal data
- pre-test: M00001 COMPLETE 10/10; M00345 2 BLOCKS MISSING 2/4

2.2. K03431 phosphoglucosamine mutase [EC:5.4.2.10]
- present in M00909
- present in gene 59755_AS26_5
    - in MAGs contigs: >59755_AS26 (MINUS_STRAND) - genomic boundaries: 4685-6080
- sequence present in Prodigal data
- pre-test: M00909 COMPLETE 5/5
 
2.3. K01965 propionyl-CoA carboxylase alpha chain [EC:6.4.1.3]
- present in M00741,M00373
- present in gene 23852_AS26_6
    - in MAGs contigs: >23852_AS26 (MINUS_STRAND) - genomic boundaries: 4984-5413
- sequence present in Prodigal data
- pre-test: M00741 COMPLETE 3/3; M00373 INCOMPLETE 5/12

3. METABAT_AS15tlH2ME_193:

3.1. K04041 fructose-1,6-bisphosphatase III [EC:3.1.3.11]
- present in M00003
- present in gene 38817_AS15_15
    - in MAGs contigs: >38817_AS15 - genomic boundaries: 12235-14215
- sequence present in Prodigal data
- pre-test: M00003 COMPLETE 8/8

3.2. K03785 3-dehydroquinate dehydratase I [EC:4.2.1.10]
- present in M00022
- present in gene 155920_AS15_4
    - in MAGs contigs: >155920_AS15 - genomic boundaries: 3365-3809
- sequence present in Prodigal data
- pre-test: M00022 COMPLETE 7/7

3.3. K00965 UDPglucose--hexose-1-phosphate uridylyltransferase [EC:2.7.7.12]
- present in M00554,M00632
- present in gene 114451_AS15_1
    - in MAGs contigs: >114451_AS15 (MINUS_STRAND) - genomic boundaries: 1-892
- sequence present in Prodigal data
- pre-test: M00632 COMPLETE 4/4; M00554 COMPLETE 2/2

4. METABAT_AS07pgkLD_51:

4.1. K01712 urocanate hydratase [EC:4.2.1.49]
- present in M00045
- present in gene 21455_AS07_6
    - in MAGs contigs: >21455_AS07 (MINUS_STRAND) - genomic boundaries: 5415-7431
- sequence present in Prodigal data
- pre-test: M00045 COMPLETE 5/5

4.2. K00768 nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase [EC:2.4.2.21]
- present in M00122
- present in gene 139372_AS07_5
    - in MAGs contigs: >139372_AS07 (MINUS_STRAND) - genomic boundaries: 2116-3232
- sequence present in Prodigal data
- pre-test: M00122 COMPLETE 7/7

4.3. K03639 GTP 3',8-cyclase [EC:4.1.99.22]
- present in M00880
- present in genes 210275_AS07_1, 79762_AS07_17
    - in MAGs contigs: >210275_AS07 - genomic boundaries: 158-818; >79762_AS07 - genomic boundaries: 19483-19936
- sequence present in Prodigal data
- pre-test: M00880 COMPLETE 5/5

5. METABAT_AS15tlH2ME_52:

5.1. K01609 indole-3-glycerol phosphate synthase [EC:4.1.1.48]
- present in M00023
- present in gene 16591_AS15_2
    - in MAGs contigs: >16591_AS15 (MINUS_STRAND) - genomic boundaries: 867-1581
- sequence present in Prodigal data
- pre-test: M00023 COMPLETE 5/5

5.2. K00767 nicotinate-nucleotide pyrophosphorylase (carboxylating) [EC:2.4.2.19]
- present in M00115
- present in gene 30510_AS15_16
    - in MAGs contigs: >30510_AS15 - genomic boundaries: 16289-17117
- sequence present in Prodigal data
- pre-test: M00115 COMPLETE 5/5

5.3. K00954 pantetheine-phosphate adenylyltransferase [EC:2.7.7.3]
- present in M00120
- present in gene 74045_AS15_15
    - in MAGs contigs: >74045_AS15 - genomic boundaries: 14178-14667
- sequence present in Prodigal data
- pre-test: M00120 5/5

6. METABAT_AS07pgkLD_55:

6.1. K00648 3-oxoacyl-[acyl-carrier-protein] synthase III [EC:2.3.1.180]
- present in M00082
- present in gene 159088_AS07_11
    - in MAGs contigs: >159088_AS07 (MINUS_STRAND) - genomic boundaries: 12306-13293
- sequence present in Prodigal data
- pre-test: M00082 COMPLETE 3/3

6.2. K01480 agmatinase [EC:3.5.3.11]
- present in M00133
- present in gene 160935_AS07_24
    - in MAGs contigs: >160935_AS07 (MINUS_STRAND) - genomic boundaries: 22723-23590
- sequence present in Prodigal data
- pre-test: M00133 COMPLETE 4/4

6.3. K02115 F-type H+-transporting ATPase subunit gamma [TC:3.A.2.1]
- present in M00157
- present in gene 126566_AS07_2
    - in MAGs contigs: >126566_AS07 - genomic boundaries: 863-1748
- sequence present in Prodigal data
- pre-test: M00157 COMPLETE 1/1

7. METABAT_AS15tlH2ME_127:

7.1. K01200 pullulanase [EC:3.2.1.41]
- present in M00855
- present in gene 29502_AS15_12
    - in MAGs contigs: >29502_AS15 (MINUS_STRAND) - genomic boundaries: 13439-15554
- sequence present in Prodigal data
- pre-test: M00855 COMPLETE 3/3

7.2. K01468 imidazolonepropionase [EC:3.5.2.7]
- present in M00045
- present in gene 38617_AS15_6
    - in MAGs contigs: >38617_AS15 - genomic boundaries: 5287-6541
- sequence present in Prodigal data
- pre-test: M00045 COMPLETE 5/5

7.3. K13038 
phosphopantothenoylcysteine decarboxylase / phosphopantothenate---cysteine ligase [EC:4.1.1.36 6.3.2.5]
- present in M00120
- present in gene 47435_AS15_2
    - in MAGs contigs: >47435_AS15 (MINUS_STRAND) - genomic boundaries: 303-909
- sequence present in Prodigal data
- pre-test: M00120 COMPLETE 5/5

8. METABAT_AS24abBPME_148:

8.1. K16370 6-phosphofructokinase 2 [EC:2.7.1.11]
- present in M00001,M00345
- present in gene 14842_AS24_3
    - in MAGs contigs: >14842_AS24 (MINUS_STRAND) - genomic boundaries: 691-1642
- sequence present in Prodigal data
- pre-test: M00001 COMPLETE 10/10; M00345 2 BLOCKS MISSING 2/4

8.2. K00872 homoserine kinase [EC:2.7.1.39]
- present in M00018
- present in gene 31761_AS24_10
    - in MAGs contigs: >31761_AS24 - genomic boundaries: 8780-9737
- sequence present in Prodigal data
- pre-test: M00018 COMPLETE 5/5

8.3. K00793 riboflavin synthase [EC:2.5.1.9]
- present in M00125
- present in gene 23060_AS24_8
    - in MAGs contigs: >23060_AS24 - genomic boundaries: 6880-7537
- sequence present in Prodigal data
- pre-test: M00125 COMPLETE 9/9

9. METABAT_AS07pgkLD_225:

9.1. K00030 isocitrate dehydrogenase (NAD+) [EC:1.1.1.41]
- present in M00009,M00010
- present in gene 171203_AS07_3
    - in MAGs contigs: >171203_AS07 - genomic boundaries: 3497-4496
- sequence present in Prodigal data
- pre-test: M00010 COMPLETE 3/3; M00009 INCOMPLETE 5/8

9.2. K00640 serine O-acetyltransferase [EC:2.3.1.30]
- present in M00021
- present in gene 9180_AS07_7
    - in MAGs contigs: >9180_AS07 (MINUS_STRAND) - genomic boundaries: 7609-8593
- sequence present in Prodigal data
- pre-test: M00021 COMPLETE 2/2

9.3. K00942 guanylate kinase [EC:2.7.4.8]
- present in M00050
- present in gene 147283_AS07_14
    - in MAGs contigs: >147283_AS07 - genomic boundaries: 11000-11618
- sequence present in Prodigal data
- pre-test: M00050 COMPLETE 4/4

10. METABAT_AS05jafATM_34:

10.1. K01695 tryptophan synthase alpha chain [EC:4.2.1.20]
- present in M00023
- present in gene 67936_AS05_1
    - in MAGs contigs: >67936_AS05 (MINUS_STRAND) - genomic boundaries: 449-1226
- sequence present in Prodigal data
- pre-test: M00023 COMPLETE 5/5

10.2. K00147 glutamate-5-semialdehyde dehydrogenase [EC:1.2.1.41]
- present in M00015
- present in gene 56766_AS05_17
    - in MAGs contigs: >56766_AS05 (MINUS_STRAND) - genomic boundaries: 21051-22302
- sequence present in Prodigal data
- pre-test: M00015 COMPLETE 3/3

10.3. K02858 3,4-dihydroxy 2-butanone 4-phosphate synthase [EC:4.1.99.12]
- present in M00125
- present in gene 17703_AS05_24
    - in MAGs contigs: >17703_AS05 - genomic boundaries: 24226-25453
- sequence present in Prodigal data
- pre-test: M00125 COMPLETE 9/9

11. METABAT_AS20ysBPTH_14:

11.1. K14941 2-phospho-L-lactate/phosphoenolpyruvate guanylyltransferase [EC:2.7.7.68 2.7.7.105]
- present in M00378
- present in gene 2020_AS20_8
    - in MAGs contigs: >2020_AS20 (MINUS_STRAND) - genomic boundaries: 9501-10149
- sequence present in Prodigal data
- pre-test: M00378 COMPLETE 4/4

11.2. K16792 methanogen homoaconitase large subunit [EC:4.2.1.114]
- present in M00608
- present in gene 65619_AS20_5
    - in MAGs contigs: >65619_AS20 (MINUS_STRAND) - genomic boundaries: 3361-4573
- sequence present in Prodigal data
- pre-test: M00608 COMPLETE 3/3

11.3. K01845 glutamate-1-semialdehyde 2,1-aminomutase [EC:5.4.3.8]
- present in M00121,M00846
- present in gene 1646_AS20_24
    - in MAGs contigs: >1646_AS20 (MINUS_STRAND) - genomic boundaries: 24749-25997
- sequence present in Prodigal data
- pre-test: M00846 COMPLETE 7/7; M00121 INCOMPLETE 7/10

12. METABAT_AS20ysBPTH_159:

12.1. K01687 dihydroxy-acid dehydratase [EC:4.2.1.9]
- present in M00019,M00570
- present in gene 45987_AS20_5
    - in MAGs contigs: >45987_AS20 (MINUS_STRAND) - genomic boundaries: 3313-4957
- sequence present in Prodigal data
- pre-test: M00019 COMPLETE 4/4; M00570 1 BLOCK MISSING 4/5

12.2. K00193 acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit beta [EC:2.3.1.169]
- present in M00357,M00422
- present in gene 72822_AS20_5
    - in MAGs contigs: >72822_AS20 - genomic boundaries: 4859-6881
- sequence present in Prodigal data
- pre-test: M00357 COMPLETE 6/6; M00422 1 BLOCK MISSING 1/2

12.3. K18209 fumarate reductase (CoM/CoB) subunit A [EC:1.3.4.1]
- present in M00620
- present in gene 41781_AS20_5
    - in MAGs contigs: >41781_AS20 - genomic boundaries: 3611-5246
- sequence present in Prodigal data
- pre-test: M00620 COMPLETE 7/7

##### Multiple ORF removal from MAG ANNOTATIONS:

Removal of the same genes indicated above from eggnog-mapper annotations (chosen also because they were identified in a single gene);
in the report_tsv (from kmc.py) the KOs of the selected genes was changed from present to absent;
in the report_txt (from kmc.py) the KOs of the selected genes were added among the missing genes.

The genes were removed using functions (rewrite_MAG_removing_genes() and others) contained in the "MAG_genomic_contig-context.ipynb" notebook. This code was also used to identify the boundaries indicated above.

The results of HMM-hits summaries (from kmc-hmm.py) were scrutinized in order to compare if the genes removed from annotations were predicted as present in the MAGs (both original and mock-MAGs, lacking part of the coding sequenes).

MAGs were grouped together according to the status of the coding sequences under test:

1. The "control" group (no removal was performed on those MAGs' contigs) is is referred to with the ID number "0".
2. The "big deletion" (70% ORF lenght) group is referred to with the ID number "1".
3. The "small deletion" (30% ORF lenght) group is referred to with the ID number "2".

#### Multiple ORF removal from Complete Genomes:

##### ORF removal (70% and 30%, respectively bigg and small deletion) in Complete Modules

1. NC_000913.3 (E. coli K-12 substr. MG1655)

1.1. K00147 glutamate-5-semialdehyde dehydrogenase [EC:1.2.1.41]
- present in M00015
- present in gene NC_000913.3_236
    - in MAGs contigs: >NC_000913.3 - genomic boundaries: 261502-262756
- sequence present in Prodigal data
- pre-test: M00015 COMPLETE 3/3

1.2. K01952 phosphoribosylformylglycinamidine synthase [EC:6.3.5.3]
- present in M00048
- present in gene NC_000913.3_2524
    - in MAGs contigs: >NC_000913.3 (MINUS_STRAND) - genomic boundaries: 2691655-2695543
- sequence present in Prodigal data
- pre-test: M00048 COMPLETE 11/11

1.3. K01637 isocitrate lyase [EC:4.1.3.1]
- present in M00012
- present in gene NC_000913.3_3927
    - in MAGs contigs: >NC_000913.3 - genomic boundaries: 4217108-4218413
- sequence present in Prodigal data
- pre-test: M00012 COMPLETE 5/5

2. NC_000908.2 (M. genitalium G37)

2.1. K00948 ribose-phosphate pyrophosphokinase [EC:2.7.6.1]
- present in M00005
- present in gene NC_000908.2_125
    - in MAGs contigs: >NC_000908.2 (MINUS_STRAND) - genomic boundaries: 66226-67318
- sequence present in Prodigal data
- pre-test: M00005 COMPLETE 1/1

2.2. K00925 phosphate acetyltransferase [EC:2.3.1.8]
- present in M00579,M00357
- present in genes NC_000908.2_778,NC_000908.2_779
    - in MAGs contigs: >NC_000908.2 (MINUS_STRAND) - genomic boundaries: 454765-455278, >NC_000908.2 (MINUS_STRAND) - genomic boundaries: 455449-455947
- sequence present in Prodigal data
- pre-test: M00048 COMPLETE 2/2

2.3. K15633 2,3-bisphosphoglycerate-independent phosphoglycerate mutase [EC:5.4.2.12]
- present in M00002,M00001,M00003
- present in gene NC_000908.2_921
    - in MAGs contigs: >NC_000908.2 (MINUS_STRAND) - genomic boundaries: 536040-537564
- sequence present in Prodigal data
- pre-test: M00002 COMPLETE 6/6; M00001 1 BLOCK MISSING 9/10 ; M00003 2 BLOCKS MISSING 6/8

3. NC_000964.3 (B. subtilis subsp. subtilis str. 168)

3.1. K00674 2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-succinyltransferase [EC:2.3.1.117]
- present in M00016
- present in gene NC_000964.3_1466
    - in MAGs contigs: >NC_000964.3 - genomic boundaries: 1488972-1489683
- sequence present in Prodigal data
- pre-test: M00016 COMPLETE 9/9

3.2. K00954 pantetheine-phosphate adenylyltransferase [EC:2.7.7.3]
- present in M00120
- present in gene NC_000964.3_1551
    - in MAGs contigs: >NC_000964.3 - genomic boundaries: 1570077-1570563
- sequence present in Prodigal data
- pre-test: M00120 COMPLETE 5/5

3.3. K08094 6-phospho-3-hexuloisomerase [EC:5.3.1.27]
- present in M00345,M00580
- present in gene NC_000964.3_351
    - in MAGs contigs: >NC_000964.3 (MINUS_STRAND) - genomic boundaries: 374602-375160
- sequence present in Prodigal data
- pre-test: M00345 COMPLETE 4/4; M00580 1 BLOCK MISSING 2/3

4. NC_009012.1 (H. thermocellum ATCC 27405)

4.1. K03786 3-dehydroquinate dehydratase II [EC:4.2.1.10]
- present in M00022
- present in gene NC_009012.1_881
    - in MAGs contigs: >NC_009012.1 (MINUS_STRAND) - genomic boundaries: 1025110-1025539
- sequence present in Prodigal data
- pre-test: M00022 COMPLETE 7/7

4.2. K00651 homoserine O-succinyltransferase/O-acetyltransferase [EC:2.3.1.46 2.3.1.31]
- present in M00017
- present in gene NC_009012.1_1924
    - in MAGs contigs: >NC_009012.1 (MINUS_STRAND) - genomic boundaries: 2188979-2189897
- sequence present in Prodigal data
- pre-test: M00017 COMPLETE 7/7

4.3. K01935 dethiobiotin synthetase [EC:6.3.3.3]
- present in M00123,M00573,M00577
- present in gene NC_009012.1_20
    - in MAGs contigs: >NC_009012.1 - genomic boundaries: 27755-28490
- sequence present in Prodigal data
- pre-test: M00123 COMPLETE 4/4; M00573 2 BLOCKS MISSING 3/5; M00577 COMPLETE 5/5

5. NZ_CP053946.1 (S. oneidensis MR-1)

5.1. K01911 o-succinylbenzoate---CoA ligase [EC:6.2.1.26]
- present in M00116
- present in gene NZ_CP053946.1_342
    - in MAGs contigs: >NZ_CP053946.1 - genomic boundaries: 390867-392313
- sequence present in Prodigal data
- pre-test: M00116 COMPLETE 9/9

5.2. K01952 phosphoribosylformylglycinamidine synthase [EC:6.3.5.3]
- present in M00048
- present in gene NZ_CP053946.1_3603
    - in MAGs contigs: >NZ_CP053946.1 (MINUS_STRAND) - genomic boundaries: 4021392-4025274
- sequence present in Prodigal data
- pre-test: M00048 COMPLETE 11/11

5.3. K02492 glutamyl-tRNA reductase [EC:1.2.1.70]
- present in M00121,M00846
- present in gene NZ_CP053946.1_4098
    - in MAGs contigs: >NZ_CP053946.1 (MINUS_STRAND) - genomic boundaries: 4574659-4575910
- sequence present in Prodigal data
- pre-test: M00121 COMPLETE 10/10; M00846 COMPLETE 7/7

##### Multiple ORF removal from Complete Genomes ANNOTATIONS:

Removal of the same genes indicated above from eggnog-mapper annotations (chosen also because they were identified in a single gene);
in the report_tsv (from kmc.py) the KOs of the selected genes was changed from present to absent;
in the report_txt (from kmc.py) the KOs of the selected genes were added among the missing genes.

The genes were removed using functions (rewrite_Genome_removing_genes() and others) contained in the "MAG_genomic_contig-context.ipynb" notebook. This code was also used to identify the boundaries indicated above.

The results of HMM-hits summaries (from kmc-hmm.py) were scrutinized in order to compare if the genes removed from annotations were predicted as present in the MAGs (both original and mock-genomes, lacking part of the coding sequenes).

Complete genomes were grouped together according to the status of the coding sequences under test:

1. The "control" group (no removal was performed on those complete genomes' contigs) is is referred to with the ID number "4".
2. The "big deletion" (70% ORF lenght) group is referred to with the ID number "5".
3. The "small deletion" (30% ORF lenght) group is referred to with the ID number "6".
