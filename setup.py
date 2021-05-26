#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import os
import re
import argparse

###############
# directories #
###############
dir_base = os.getcwd()+"/"

Modules_directory = dir_base+"KEGG_MODULES_ref/"
kkfiles_directory = Modules_directory+"file_test_kk/"
KAnnotation_directory = dir_base+"KEGG_mappings/"
ktests_directory = dir_base+"ktests/"
klists_directory = dir_base+"klists/"
report_directory = dir_base+"reports/"
report_tsv_directory = dir_base+"reports_tsv/"

def set_directories(dir_base, hmm_usage = False, gapfill_usage = False):
    directories_to_make = []
    os.chdir(dir_base)
    initial_list = os.listdir()

    KAnnotation_directory = "KEGG_mappings"
    ktests_directory = "ktests"
    klists_directory = "klists"
    report_txt_directory = "reports_txt"
    report_tsv_directory = "reports_tsv"

    directories_to_make.append(KAnnotation_directory)
    directories_to_make.append(ktests_directory)
    directories_to_make.append(klists_directory)
    directories_to_make.append(report_txt_directory)
    directories_to_make.append(report_tsv_directory)

    if hmm_usage:
        taxa_dir = "taxonomies"
        dir_base_KO = "Knumber_ntsequences"
        dir_base_KOaa = "Knumber_aasequences"
        msa_dir = "multiple_fasta"
        hmm_dir = "HMM"
        hmm_hits_dir = "HMM_HITS"
        dir_genomes = "genomes"
        oneBM_modules_dir = "oneBM_modules"

        directories_to_make.append(taxa_dir)
        directories_to_make.append(dir_base_KO)
        directories_to_make.append(dir_base_KOaa)
        directories_to_make.append(msa_dir)
        directories_to_make.append(hmm_dir)
        directories_to_make.append(hmm_hits_dir)
        directories_to_make.append(dir_genomes)
        directories_to_make.append(oneBM_modules_dir)

        os.system("curl --silent http://rest.kegg.jp/get/br:br08601 > br08601.keg")
        print("KEGG Organisms hierarchy DOWNLOADED")

        if not "genomes.instruction" in os.listdir():
            os.system("echo 'id\ttaxonomy\tuniverse\n' > genomes.instruction")
            print("genome_instruction file GENERATED")
        else:
            print("Instruction file already EXISTS")
        if not "module_file.instruction" in os.listdir():
            os.system("touch module_file.instruction")
            print("module_file GENERATED")
        else:
            print("module_file already EXISTS")
        if not "ko_file.instruction" in os.listdir():
            os.system("touch ko_file.instruction")
            print("ko_file GENERATED")
        else:
            print("ko_file already EXISTS")

    if gapfill_usage:
        gapfill_report_directory = "report_gapfill"
        bigg_api = "biggapi_download"
        DB_directory = "DB"
        model_directory = "models"
        gapfilled_model_directory = "models_gapfilled"
        de_novo_model_directory = "de_novo_models"

        directories_to_make.append(gapfill_report_directory)
        directories_to_make.append(bigg_api)
        directories_to_make.append(DB_directory)
        directories_to_make.append(model_directory)
        directories_to_make.append(gapfilled_model_directory)
        directories_to_make.append(de_novo_model_directory)

    for el in directories_to_make:
        if not el in initial_list:
            os.system("mkdir "+el)
            print(el+" folder CREATED")
        elif el in initial_list:
            print(el+" already EXISTS")

def set_kk_database():
    pass

def update_kk_database():
    pass

###############################################################################

def main():
    parser = argparse.ArgumentParser(description='''Setup command for KEMET analysis.
                                                    Create folders and generate/update KEGG Module .kk database''')

    parser.add_argument('-k','--set_kk_DB', action = "store_true",
                        help='''Choose this option in order to create KEGG Module DB (.kk files), in order to use KMC.''') #TODO
    parser.add_argument('-u','--update_kk_DB', action = "store_true",
                        help='''Choose this option in order to update already existing KEGG Module DB (.kk files).''') #TODO

    parser.add_argument('-H','--hmm_usage', action = "store_true",
                        help='''Choose this option in order to create and download required folders for HMM analysis, follow-up of KMC.''')
    parser.add_argument('-G','--gapfill_usage', action = "store_true",
                        help='''Choose this option in order to create and download required folders for Gapfilling, follow-up of KMC+HMM.''')
    parser.add_argument('-v','--verbose', action ="store_true",
                        help='''Print more informations - for debug or log.''') #TODO: add verbosity to commands
    args = parser.parse_args()

###############################################################################

    if args.hmm_usage:
        if args.gapfill_usage:
            set_directories(dir_base, hmm_usage=True, gapfill_usage=True)
        else:
            set_directories(dir_base, hmm_usage=True)
    else:
        set_directories(dir_base)
    #TODO: ADD THE DATABASE OF kk-files WITH THE UPDATED MODULES LIST
    if args.set_kk_DB:
        set_kk_database()
    if args.update_kk_DB:
        update_kk_database()

if __name__ == "__main__":
    main()
