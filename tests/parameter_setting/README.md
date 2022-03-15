# Setting a parametric threshold for quality hits recovery  

Ideally, setting a threshold value had the purpose of minimizing low quality HMM hits e.g. genomic sequences corresponding to only a portion of the protein of interest (domains), or hits that correspond to different KEGG Orthologs (KOs).  
At the same time, the predictive power needs to be useful to justify the adoption of KEMET as a way to enhance functional annotation, i.e. maximixing high quality HMM hits.  

HMMs performances were first leveraged using `kemet.py` without any score threshold indication (`--threshold_value 0.0`) for different KOs subgroup in two distinct MAGs dataset.  
These represented two subsets of high-quality MAGs from two recently published papers (10.1016/j.scitotenv.2021.146296 and 10.1186/s40168-019-0780-9).  
Specifically, the two KO subgroups were:  
- KOs missing from “Carbohydrate metabolism” Modules annotations in the first MAG dataset  

- all KOs missing from 1 block missing Modules in the second MAG dataset  

Sequences recovered with KEMET were then manually inspected using BLASTp against NCBI nr database, assigning a quality value to HMM hits.  

The now default `kemet.py` threshold selected by the Authors (0.43), had the following for the two test dataset:  

- KEMET identifies 48/52 and 124/131 validated orthologs for the two datasets, which are not included in the original KO functional annotation (assigned with eggnog-mapper v2, in March 2021.  

- Conversely, 4/52 and 7/131 KO assignments via HMM result in bad quality hits (i.e. hits resulting in different functional identification from BLASTp), and thus return 7.5% and 4.6% false positive rates, respectively.  Tabular results of these KEMET analyses are included (`KEMET/tests/raw/`) alongside hits quality in the present folder.