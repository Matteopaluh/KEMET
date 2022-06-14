#!/usr/bin/env python
# coding: utf-8

import os, re
from os import path
import argparse

_fasta_extensions = [".fa", ".fna", ".fasta"]

def kegg_taxonomy_from_gtdbtk_NCBI_mapping(ncbi_from_gtdb_file, genomes_info_file, fasta_extension, verbose=False):
    rank_prefixes = ["d__","p__", "c__", "o__", "f__", "g__", "s__"]
    excluded_phyla = [
        "Abditibacteriota",
        "Armatimonadetes",
        "Atribacterota",
        "Balneolaeota",
        "Caldiserica",
        "Calditrichaeota",
        "Chrysiogenetes",
        "Coprothermobacterota",
        "Dictyoglomi",
        "Fibrobacteres",
        "Ignavibacteriae",
        "Kiritimatiellaeota",
        "Lentisphaerae",
        "Nanoarchaeota",
        "Nitrospinae",
        "Rhodothermaeota",
    ]

    ncbi_kegg_common_phyla = [
        "Acidobacteria",
        "Actinobacteria",
        "Aquificae",
        "Bacteroidetes",
        "Chlamydiae",
        "Chlorobi",
        "Chloroflexi",
        "Crenarchaeota",
        "Cyanobacteria",
        "Deferribacteres",
        "Deinococcus-Thermus",
        "Elusimicrobia",
        "Euryarchaeota",
        "Fusobacteria",
        "Gemmatimonadetes",
        "Nitrospirae",
        "Planctomycetes",
        "Spirochaetes",
        "Synergistetes",
        "Tenericutes",
        "Thaumarchaeota",
        "Thermodesulfobacteria",
        "Thermotogae",
        "Verrucomicrobia",
    ]

    ncbi_to_kegg_mapping = {}
    genomes_lacking_info = []
    genomes_lacking_kegg_representatives = []

    with open(ncbi_from_gtdb_file) as f:
        for linum, line in enumerate(f):
            line = line.strip().split("\t")
            if linum == 0:
                ncbi_header = line.index("NCBI classification")
                continue

            genome = line[0]
            ncbi_taxonomy = line[ncbi_header]
            for rank in rank_prefixes:
                ncbi_taxonomy = re.sub(rank, "", ncbi_taxonomy)
            _domain, _phylum, _class, _order, _family, _genus, _species = ncbi_taxonomy.split(";")
            
            if _domain == "" or _phylum == "":
                genomes_lacking_info.append(genome)
                continue

            if _phylum in excluded_phyla:
                genomes_lacking_kegg_representatives.append(genome)
                continue

            if _phylum in ncbi_kegg_common_phyla:
                ncbi_to_kegg_mapping.update({genome : _phylum})

            if _phylum == "Proteobacteria":
                if (_class == "Alphaproteobacteria"
                    or _class == "Betaproteobacteria"
                    or _class == "Deltaproteobacteria"
                    or _class == "Epsilonproteobacteria"
                    ):
                    ncbi_to_kegg_mapping.update({genome : _class})
                elif _class == "Gammaproteobacteria":
                    if _order == "Enterobacteriales":
                        ncbi_to_kegg_mapping.update({genome : "Gammaproteobacteria - Enterobacteria"})
                    else:
                        ncbi_to_kegg_mapping.update({genome : "Gammaproteobacteria - Others"})
                else:
                    ncbi_to_kegg_mapping.update({genome : "Other Proteobacteria"})

            if _phylum == "Firmicutes":
                if (_class == "Bacilli"
                    or _class == "Clostridia"
                    ):
                    ncbi_to_kegg_mapping.update({genome : f"{_phylum} - {_class}"})
                else:
                    ncbi_to_kegg_mapping.update({genome : "Firmicutes - Others"})
    if verbose:
        for genome in genomes_lacking_info:
            print(f"{genome} not included: minimum taxonomy information lacking. Check script help page (-h).")

        for genome in genomes_lacking_kegg_representatives:
            print(f"{genome} not included: not enough KEGG Organisms with the same phylum. Check script help page (-h).")

    if path.isfile(genomes_info_file):
        with open(genomes_info_file, "a") as g:
            for genome, kegg_taxonomy in ncbi_to_kegg_mapping.items():
                print(genome+fasta_extension, kegg_taxonomy, "", file=g, sep="\t")

        print("The {} file has been updated with {} genome(s) taxonomy indications, using '{}' extension.".format(genomes_info_file, str(len(ncbi_to_kegg_mapping)), fasta_extension))      
    else:
        print(f"""{genomes_info_file} path does not lead to KEMET 'genomes.instruction' file.
            Be sure to use the proper command to set KEMET working directory first.
            If this message still gets printed, double check the indicated path.
            """)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description=
    '''
    Add necessary taxonomy informations of MAGs/Genomes of interest for KEMET HMM and GSMM analyses.
    Use this after the GTDB-tk "gtdb_to_ncbi_majority_vote.py" script (that converts from GTDB taxonomy to NCBI),
    to further convert to KEGG BRITE taxonomy.
    This script will include info on the MAGs/Genomes indicated in the output file from the aforementioned script.
    IMPORTANT:

    The automatic taxonomy conversion has notable exceptions, such as "Candidate" phyla,
    as well as other phyla lacking sufficient (3+) KEGG Organism representatives.
    ''')
    parser.add_argument('-i','--add_genomes_instruction_file', required = True,
                        help='''Include the relative path to KEMET "genomes.instruction" file.''')
    parser.add_argument('-t','--add_gtdb_to_ncbi_output', required = True,
                        help='''Include the relative path to "gtdb_to_ncbi_majority_vote.py" output file.''')
    parser.add_argument('-f','--fasta_extension', required = True, choices=_fasta_extensions,
                        help='''Complete "genomes.instruction" file names with the indicated extension.''')
    parser.add_argument('-v','--verbose', action ="store_true",
                        help='''Print more informations - for debug or log.''')
    args = parser.parse_args()

###############################################################################

    genomes_info_file = args.add_genomes_instruction_file
    ncbi_from_gtdb_file = args.add_gtdb_to_ncbi_output
    fasta_extension = args.fasta_extension

    if args.verbose:
        kegg_taxonomy_from_gtdbtk_NCBI_mapping(ncbi_from_gtdb_file, genomes_info_file, fasta_extension, verbose=True)
    else:
        kegg_taxonomy_from_gtdbtk_NCBI_mapping(ncbi_from_gtdb_file, genomes_info_file, fasta_extension)

if __name__ == "__main__":
    main()