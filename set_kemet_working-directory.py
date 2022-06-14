#!/usr/bin/env python
# coding: utf-8

import os
from os import path
import argparse
#import shutil
#from kemet_data import project_dir

def set_directories(dir_base, gapfill_usage=False):
    """
    Setup function to generate folders and instruction files for
    KEMET script execution.

    Args:
        dir_base                   (str): folder path in which "kemet.py" is going to be executed
        gapfill_usage   (bool, optional): flag to indicate that GSMMs functionalities are wanted.
                                          Defaults to False.
    """

    directories_to_make = [
        "KEGG_annotations",
        "ktests",
        "klists",
        "reports_txt",
        "reports_tsv",
        "taxonomies",
        "Knumber_ntsequences",
        "multiple_fasta",
        "HMM",
        "HMM_HITS",
        "genomes",
        "oneBM_modules",
    ]

    if gapfill_usage:
        directories_to_make += [
            "report_gapfill",
            "biggapi_download",
            "models",
            "models_gapfilled",
            "de_novo_models",
        ]

    genome_instruction_file = path.join(dir_base, "genomes.instruction")
    module_file = path.join(dir_base, "module_file.instruction")
    ko_file = path.join(dir_base, "ko_file.instruction")
    kegg_brite_organisms = path.join(dir_base, "br08601.keg")

    os.system(f"curl --silent https://rest.kegg.jp/get/br:br08601 > {kegg_brite_organisms}")
    print("KEGG Organisms hierarchy DOWNLOADED")

    if path.isfile(genome_instruction_file):
        print("Instruction file ALREADY EXISTS")
    else:
        with open(genome_instruction_file, "w") as f:
            print("id", "taxonomy", "universe", sep="\t", file=f)
        print("genome_instruction file GENERATED")
    
    if path.isfile(module_file):
        print("module_file ALREADY EXISTS")
    else:
        os.mknod(module_file)
        print("module_file GENERATED")
        
    if path.isfile(ko_file):
        print("ko_file ALREADY EXISTS")
    else:
        os.mknod(ko_file)
        print("ko_file GENERATED")

    for el in directories_to_make:
        dir = path.join(dir_base, el)
        if path.isdir(dir):
            print(f"{dir} folder ALREADY EXISTS")
        else:
            os.mkdir(dir)
            print(f"{dir} folder CREATED")

#def copy_ref_kegg_modules_and_DB(project_dir, dir_base):
#    if not path.isdir("kemet_data"):
#        os.mkdir("kemet_data")
#        os.chdir("kemet_data")
#        shutil.copytree(project_dir, path.join(dir_base, "kemet_data"))

def set_kk_database():
    pass

def update_kk_database():
    pass

###############################################################################

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description=
    '''
    Base command for setting KEMET package working directory.
    Create folders and instruction files; helper function to manage KEGG MODULE .kk database.
    ''')

    # Add command line option for setting the base directory?
    parser.add_argument('-k','--set_kk_DB', action="store_true",
                        help='''
                        Choose this option to generate KEGG Module DB (.kk files),
                        in order to perform KEGG Modules Completeness evaluation.
                        Default: already generated''')
    parser.add_argument('-u','--update_kk_DB', action="store_true",
                        help='''
                        Choose this option to update already existing KEGG Module DB (.kk files).''')
    parser.add_argument('-G','--gapfill_usage', action="store_true",
                        help='''
                        Choose this option to create required folders for the GSMM Gapfilling,
                        follow-up of the HMM search procedures.''')
    args = parser.parse_args()

###############################################################################

    dir_base = os.getcwd()

    set_directories(dir_base, args.gapfill_usage)

    #POSSIBILITY: ALLOW DYNAMIC kk-files DATABASE UPDATES
    if args.set_kk_DB:
        set_kk_database()

    elif args.update_kk_DB:
        update_kk_database()

    #if not "kemet_data" in os.listdir() and not "KEGG_MODULES" in os.listdir():
    #    copy_ref_kegg_modules_and_DB(project_dir, dir_base)

if __name__ == "__main__":
    main()