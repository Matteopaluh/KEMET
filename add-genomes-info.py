#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import re
import argparse

def add_organism_instruction(org_id, taxonomy, universe, verbose=False):
    line = org_id+"\t"+taxonomy+"\t"+universe+"\n"
    if not "genomes.instruction" in os.listdir():
        print("ERROR: genomes.instruction file missing - RUN setup.py")
        return 0
    else:
        f = open("genomes.instruction", "a")
        f.write(line)
        f.close()
        print("COMPLETE organism instruction:\n"+line.strip().replace("\t", " "))

###############################################################################

def main():
    parser = argparse.ArgumentParser(description='''Add instructions for MAGs/Genomes of interest for KMC+HMM analysis''')

    parser.add_argument('-g','--add_fasta_file', required = True,
                        help='''Add info regarding MAG/Genome file name (without extension).''')
    parser.add_argument('-t','--add_taxonomy', required = True,
                        help='''Add info regarding MAG/Genome taxonomy.''') #TODO: PARSE THE INFO FROM GTDB TAXONOMY FILE?
    parser.add_argument('-u','--add_universe', default = "",
                        help='''Add info regarding the metabolic universe to which the organism belong (WHEN TO USE: further metabolic reconstruction).''')

    parser.add_argument('-v','--verbose', action ="store_true",
                        help='''Print more informations - for debug or log.''')
    args = parser.parse_args()

###############################################################################

    org_id = args.add_fasta_file
    taxonomy = args.add_taxonomy
    universe = args.add_universe

    if args.verbose:
        add_organism_instruction(org_id, taxonomy, universe, verbose=True)
    else:
        add_organism_instruction(org_id, taxonomy, universe)

if __name__ == "__main__":
    main()

