#!/usr/bin/env python
# coding: utf-8

# Imports
import os
import re
from multiprocessing import Pool
from datetime import datetime
import argparse
import reframed
from reframed import load_cbmodel
from reframed import save_cbmodel

###############
# extra specs #
###############
_ktest_formats = ["eggnog", "kaas", "kofamkoala"]
_hmm_modes = ["onebm","modules","kos"]
_def_thr = 0.43 # threshold optimized from test datasets
_gapfill_modes = ["existing","denovo"]
_base_com_KEGGget = "curl --silent http://rest.kegg.jp/get/"

# External dependencies base commands - experienced users can edit variables' with proper parameters e.g. to modify threads etc.
_base_com_mafft = "mafft --quiet --auto --thread -1 MSA_K_NUMBER.fna > K_NUMBER.msa"
_base_com_hmmbuild = "hmmbuild --informat afa K_NUMBER.hmm K_NUMBER.msa > /dev/null"
_base_com_nhmmer = "nhmmer --tblout K_NUMBER.hits K_NUMBER.hmm PATHFILE > /dev/null"

def _timeinfo():
    """
    Helper function for generating time indications on the console,
    to keep track of script process.

    Returns:
    timeinfo        (str): present date and time, up to seconds
    """
    timeinfo = str(datetime.now().strftime("%Y-%m-%d %H-%M-%S"))+"\t"
    return timeinfo

def eggnogXktest(eggnog_file, converted_output, KAnnotation_directory, ktests_directory):
    """
    Starting from eggNOG pre-annotations (".emapper.annotations"),
    (1 gene - many annotations), keeps only KOs in a converted ".ktest" file.

    Args:
        eggnog_file             (str): eggnog mapper output file name
        converted_output        (str): resulting intermediate .ktest file name
        KAnnotation_directory   (str): eggnog mapper output folder path
        ktests_directory        (str): output .ktest file folder path

    Returns:
        converted_output        (str): output .ktest file name
        KOs                    (dict): (not used in current version)
    """
    os.chdir(KAnnotation_directory)
    KOs = {}

    with open(eggnog_file) as g:
        for linum, line in enumerate(g, -1):
            if not line.startswith("#"):
                break
        g.seek(0)
        headers = g.readlines()[linum].strip().split("\t")
        koslice = 0
        for field in headers:
            if not field == "KEGG_ko":
                koslice += 1
            else:
                break
        g.seek(0)
        for line in g.readlines():
            if not line.startswith("#"): # skip header & info lines w/o genes
                #fasta_id = line.strip().split("\t")[0]
                egg_kos = line.strip().split("\t")[koslice].replace("ko:","")
                if egg_kos != "" and egg_kos != "-":
                    egg_kos_hits = egg_kos.split(",")
                    for ko in egg_kos_hits:
                        if not ko in KOs:
                            KOs[ko] = 1
                        else:
                            KOs[ko] += 1

            # POSSIBILITY: for each gene, correcting per diff. ortholog hits - if more KOs -> fraction of KO count
                        #if not ko in KOs:
                            #KOs[ko] = round(1/len(egg_kos_hits), 2)
                        #else:
                            #KOs[ko] += round(1/len(egg_kos_hits), 2)
                else:
                    pass

    try:
        os.chdir(ktests_directory)
    except:
        os.mkdir(ktests_directory)
        os.chdir(ktests_directory)
    g = open(converted_output, "w")
    for ko in KOs.keys():
        g.write(ko+"\n")
    g.close()

    return converted_output, KOs

def KAASXktest(file_kaas, converted_output, KAnnotation_directory, ktests_directory):
    """
    Starting from KAAS output (1 gene - 1 KO),
    keeps only KOs in a converted ".ktest" file.

    Args:
        file_kaas               (str): KAAS output file name
        converted_output        (str): resulting intermediate .ktest file name
        KAnnotation_directory   (str): KAAS output folder path
        ktests_directory        (str): output .ktest file folder path

    Returns:
        converted_output        (str): output .ktest file name
        KOs                    (dict): (not used in this version)
    """
    os.chdir(KAnnotation_directory)
    KOs = []
    with open(file_kaas) as f:
        for line in f.readlines():
            line_s = line.strip().split("\t")
            if len(line_s) == 2:
                if not line_s[1] in KOs:
                    KOs.append(line_s[1])
            elif len(line_s) == 1:
                continue
    try:
        os.chdir(ktests_directory)
    except:
        os.mkdir(ktests_directory)
        os.chdir(ktests_directory)
    g = open(converted_output, "w")
    for ko in KOs:
        g.write(ko+"\n")
    g.close()
    return converted_output, KOs

def kofamXktest(kofamkoala_file, converted_output, KAnnotation_directory, ktests_directory):
    """
    Starting from KofamKOALA output (1 gene - many annotations),
    keeps only KOs in a converted ".ktest" file.

    Args:
        kofamkoala_file         (str): KofamKOALA output file name
        converted_output        (str): resulting intermediate .ktest file name
        KAnnotation_directory   (str): KofamKOALA output folder path
        ktests_directory        (str): output .ktest file folder path

    Returns:
        converted_output        (str): output .ktest file name
        KOs                    (dict): (not used)
    """
    os.chdir(KAnnotation_directory)
    KOs = {}

    with open(kofamkoala_file) as g:
        # look for fasta_id & KO info based on gene calling IDs lenght
        spacer = g.readlines()[1].strip()
        fastaslice = spacer.index(" ",1)+1
        koslice = spacer.index(" ",fastaslice)+1
        g.seek(0)
        
        for line in g.readlines()[2:]: # skip header and spacer lines w/o genes
            fasta_id = line[:fastaslice].strip().replace("* ","")
            kofam_ko = line[fastaslice:koslice].strip()
            if kofam_ko != "":
                if not kofam_ko in KOs:
                    KOs[kofam_ko] = 1
                else:
                    KOs[kofam_ko] += 1

            else:
                pass
    try:
        os.chdir(ktests_directory)
    except:
        os.mkdir(ktests_directory)
        os.chdir(ktests_directory)
    g = open(converted_output, "w")
    for ko in KOs.keys():
        g.write(ko+"\n")
    g.close()
    return converted_output, KOs

def create_KO_list(file_ko_list, ktests_directory):
    """
    Returns a Python list-object from KOs file as recovered in pre-annotations (".ktest" file).

    Args:
        file_ko_list            (str): KO list file name (".ktest")
        ktests_directory        (str): ".ktest" input file folder path

    Returns:
        ko_list                (list): list of single KOs present in the input pre-annotation
    """
    os.chdir(ktests_directory)
    ko_list = []
    with open(file_ko_list) as f:
        for line in f.readlines():
            line_s = line.strip()
            ko_list.append(line_s)
    return ko_list

def testcompleteness(ko_list, kk_file, kkfiles_directory, report_txt_directory, file_output, cutoff = 0):
    """
    Computes KEGG Modules completeness from KO pre-annotation.
    Reports that in a flat-file,
    including single missing KOs and their position, relative to KEGG Modules blocks

    Args:
        ko_list                (list): output from previous "create_KO_list" function
        kk_file                 (str): .kk file, which includes a KEGG Module definition properly parsed and organized in blocks
        kkfiles_directory       (str): input .kk files folder
        report_txt_directory    (str): output folder
        file_output             (str): output file name
        cutoff        (int, optional): Minimum obtained KEGG Module completeness percentage to be included in report file.
                                        Defaults to 0.
    """
    os.chdir(kkfiles_directory)
    report = []
    with open(kk_file) as f:
        count_lines = 0
        linenumber = 0
        missing = 0
        present = 0
        complexes = []
        optional = []
        submodules = False
        submodule = ""
        subOR_presence = False
        subOR_dict = {}
        to_remove = []
        ko_list_optional = ko_list
        extended_name = f.readline().strip().replace(".txt","")
        f.seek(0)
        report.append(kk_file+"\t"+extended_name+"\n")
        v = f.readlines()
        end = len(v)

# search presence-absence complexes//optional//list-type
        f.seek(0)
        if "COMPLEXES_LIST\n" in v:
            v_complex = v.index("COMPLEXES_LIST\n")
        else:
            v_complex = -1
        if "OPTIONAL_LIST\n" in v:
            v_optional = v.index("OPTIONAL_LIST\n")
        else:
            v_optional = -1
        if "/\n" in v:
            submodules = True
            submodule = 0
        if "//\n" in v:
            submodules = True
            submodule = 0
            subOR_presence = True
            subOR = -1

# search complexes
        if v_complex != -1:
            end = v_complex
            c_list = v[v_complex+1].replace("\n", "").replace("\t", "").split(", ")
            for el in c_list:
                complexes.append(el)
# search optionals
        if v_optional != -1:
            o_list = v[v_optional+1].replace("\n", "").replace("\t", "").split(", ")
            for el in o_list:
                optional.append(el)
            ko_list_optional = [el for el in ko_list]
            for el in optional:
                ko_list_optional.append(el)

        for line in v[1:end]:
            count_lines += 1
            if submodules:
                if line == "/\n":
                    linenumber = 0
                    check = 0
                    submodule += 1
                    continue
                if subOR_presence:
                    end_or = end-1
                    if line == "//\n"  or count_lines == end_or:
                        linenumber = 0
                        check = 0
                        submodule += 1
                        if not count_lines == end_or:
                            subOR += 1
                        if subOR != 0:
                            if not count_lines == end_or:
                                tmp = str(present)+"__"+str(present+missing)
                                subOR_infos += tmp
                                subOR_dict.update({subOR:subOR_infos})
                                present = 0
                                missing = 0

                        subOR_infos = ""
                        if not count_lines == end_or:
                            continue
                        else:
                            pass

            linenumber += 1
            check = 0
            ko_line = line.strip().split(", ")

            if len(complexes) != 0:
                while check == 0:
                    for singlecomplex in complexes:
                        k_singlecomplex = re.split("[+-]", singlecomplex.strip())
                        if all(el in ko_line for el in k_singlecomplex):
                            if all(el in ko_list_optional for el in k_singlecomplex):
                                check = 1
                            else:
                                continue
                    else:
                        for element in ko_line:
                            if not element in str(complexes):
                                if element in ko_list:
                                    check = 1
                                    break
                                else:
                                    continue
                        else:
                            break

                if check == 1:
                    present += 1
                elif check == 0:
                    missing += 1
                    control = str(linenumber)+"."+str(submodule)+"\t"+str(line)
                    report.append(control)
                    pass

            elif len(complexes) == 0:
                for element in ko_list:
                    if element in line:
                        check = 1
                else:
                    if check == 1:
                        present += 1
                        pass
                    else:
                        missing += 1
                        control = str(linenumber)+"."+str(submodule)+"\t"+str(line)
                        report.append(control)
                        pass

            total = present+missing
            percentage_round = round((present/(total))*100, 2)

        else:
            if subOR_presence and count_lines == end_or:
                subOR += 1
                tmp = str(present)+"__"+str(present+missing)
                subOR_infos += tmp
                subOR_dict.update({subOR:subOR_infos})
                percentage_round = -1

        for subOR, value in subOR_dict.items():
            tmp_present = int(value.split("__")[0])
            tmp_total = int(value.split("__")[1])
            tmp_percentage_round = round((tmp_present/(tmp_total))*100,2)

            if tmp_percentage_round > percentage_round:
                present = tmp_present
                total = tmp_total
                percentage_round = round((present/(total))*100,2)
                subOR_most = subOR
        
        if percentage_round == 100:
            completeness = "COMPLETE"
        else:
            completeness = "INCOMPLETE"
        report.insert(1, "%\t"+str(percentage_round)+"\t"+str(present)+"__"+str(total)+"\t"+completeness+"\n")

        if subOR_presence:
            for info in report[2:]:
                info = info.strip()
                sub = int(info.split(".")[1].split("\t")[0].strip())
                if sub != subOR_most:
                    to_remove.append(info)

        for el in to_remove:
            report.remove(el+"\n")

    if percentage_round >= cutoff:
        try:
            os.chdir(report_txt_directory)
        except:
            os.mkdir(report_txt_directory)
            os.chdir(report_txt_directory)
        g = open(file_output, "a")
        for el in report:
            g.write(el)
        g.write("\n")
        g.close()

def testcompleteness_tsv(ko_list, kk_file, kkfiles_directory, report_tsv_directory, file_report_tsv, as_kegg = False, cutoff = 0):
    """
    Computes KEGG Modules completeness from KO pre-annotation.
    Reports that in a tab-separated file,
    including single missing KOs and their position, relative to KEGG Modules blocks.

    Args:
        ko_list                    (list): output from previous "create_KO_list" function
        kk_file                     (str): .kk file, which includes a KEGG Module definition properly parsed and organized in blocks
        kkfiles_directory           (str): input .kk files folder
        report_tsv_directory        (str): output folder path
        file_report_tsv             (str): output file name
        as_kegg                    (bool): option to report KEGG Modules completeness as KEGG mapper (see README for details)
        cutoff            (int, optional): Minimum obtained KEGG Module completeness percentage to be included in report file.
                                            Defaults to 0.
    """
    os.chdir(kkfiles_directory)
    report = []
    report_tsv = []
    with open(kk_file) as f:
        linenumber = 0
        count_lines = 0
        missing = 0
        present = 0
        KOmodule = []
        Kmissing = []
        Kpresent = []
        complexes = []
        optional = []
        submodules = False
        subAND_presence = False
        subOR_presence = False
        subOR_dict = {}
        ko_list_optional = ko_list
        extended_name = f.readline().strip().replace(".txt","")
        f.seek(0)
        report.append(kk_file+"\t"+extended_name+"\t")
        report_tsv.append(kk_file[:-3])
        report_tsv.append(extended_name[7:])
        v = f.readlines()
        end = len(v)
        f.seek(0)
    # flags complexes/optionals/list-indications presence or absence

        #COMPLEXES
        if "COMPLEXES_LIST\n" in v:
            v_complex = v.index("COMPLEXES_LIST\n")
        else:
            v_complex = -1

        #OPTIONALS
        if "OPTIONAL_LIST\n" in v:
            v_optional = v.index("OPTIONAL_LIST\n")
        else:
            v_optional = -1

        #LIST-INDICATIONS
        if "/\n" in v:
            submodules = True
            subAND_presence = True
            submodule = 0
            subAND = 0
        if "//\n" in v:
            submodules = True
            subOR_presence = True
            submodule = 0
            subOR = -1

    # list complexes
        if v_complex != -1:
            #KOs search stops at complex line
            end = v_complex
            c_list = v[v_complex+1].replace("\n", "").replace("\t", "").split(", ")
            for el in c_list:
                complexes.append(el)

    # list optionals
        if v_optional != -1:
            o_list = v[v_optional+1].replace("\n", "").replace("\t", "").split(", ")
            for el in o_list:
                optional.append(el)
            ko_list_optional = [el for el in ko_list]
            for el in optional:
                ko_list_optional.append(el)

    # KOs presence/absence
        for line in v[1:end]:
            ko_line = line.strip().split(", ")
            for single_ko in ko_line:
                if single_ko == "/" or single_ko == "//":
                    continue
                KOmodule.append(single_ko)

        for KO in KOmodule:
            if KO in ko_list:
                if not KO in Kpresent:
                    Kpresent.append(KO)
            else:
                Kmissing.append(KO)

    # CHECKS for each line in .kk file: KOs, complexes, optionals and list-indication
        for line in v[1:end]:
            count_lines += 1
            check = 0
            linenumber += 1

            if submodules:
                if line == "/\n":
                    submodule += 1
                    subAND += 1
                    linenumber = 0
                    continue
                if subOR_presence:
                    end_or = end-1
                    if line == "//\n" or count_lines == end_or:
                        submodule += 1
                        subOR += 1
                        linenumber = 0 
                        if subOR != 0:
                            if not count_lines == end_or:
                                tmp = str(present)+"__"+str(present+missing)
                                subOR_infos += tmp
                                subOR_dict.update({subOR:subOR_infos})
                                present = 0
                                missing = 0

                        subOR_infos = ""
                        if not count_lines == end_or:
                            continue
                        else:
                            pass

            ko_line = line.strip().split(", ")
            if len(complexes) != 0:
                while check == 0:
                    for singlecomplex in complexes:
                        k_singlecomplex = re.split("[+-]", singlecomplex.strip())
                        if all(el in ko_line for el in k_singlecomplex):
                            if all(el in ko_list_optional for el in k_singlecomplex):
                            # this way: if EACH KO of complex is present in genome KOs, CHECK positive!
                                check = 1
                            else:
                                continue
                    else:
                        for ko in ko_line:
                            if not ko in str(complexes):
                                if ko in ko_list:
                                    check = 1
                                    break
                                else:
                                    continue
                        else:
                            break

                if check == 1:
                    present += 1
                elif check == 0:
                    missing += 1
                    pass

            elif len(complexes) == 0:
                for ko in ko_list:
                    if ko in line:
                        check = 1
                else:
                    if check == 1:
                        present += 1
                        pass
                    else:
                        missing += 1
                        pass

            total = present+missing
            missing_blocks = str(present)+"__"+str(total)
            percentage_round_tsv = round((present/(total))*100,2)
        else:
            if subOR_presence and count_lines == end_or:
                tmp = str(present)+"__"+str(present+missing)
                subOR_infos += tmp
                subOR_dict.update({subOR:subOR_infos})

    # check better completeness from alternative sub-modules
        if subOR_presence:
            for value in subOR_dict.values():
                tmp_present = int(value.split("__")[0])
                tmp_total = int(value.split("__")[1])
                tmp_percentage_round_tsv = round((tmp_present/(tmp_total))*100,2)

                if tmp_percentage_round_tsv > percentage_round_tsv:
                    present = tmp_present
                    total = tmp_total
                    missing_blocks=str(present)+"__"+str(total)
                    percentage_round_tsv = round((present/(total))*100,2)

        if not as_kegg:
            if present == total:
                completeness_tsv = "COMPLETE"
            if present+2 < total:
                completeness_tsv = "INCOMPLETE"
            elif present+2 == total:
                completeness_tsv = "2 BLOCKS MISSING"
            elif present+1 == total:
                completeness_tsv = "1 BLOCK MISSING"
        if as_kegg:
            if present == total:
                completeness_tsv = "COMPLETE"
            elif present+2 < total or total < 3:
                completeness_tsv = "INCOMPLETE"
            elif present+2 == total:
                completeness_tsv = "2 BLOCKS MISSING"
            elif present+1 == total:
                completeness_tsv = "1 BLOCK MISSING"

        report_tsv.append(completeness_tsv)
        report_tsv.append(missing_blocks)
        report_tsv.append(Kmissing)
        report_tsv.append(Kpresent)

    # IO-files operations
    try:
        os.chdir(report_tsv_directory)
    except:
        os.mkdir(report_tsv_directory)
        os.chdir(report_tsv_directory)

    if percentage_round_tsv >= cutoff: # OPTIONAL
        h = open(file_report_tsv, "a")
        for el in report_tsv:
            if type(el) == str:
                h.write(el+"\t")
            if type(el) == list:
                knums = ",".join(el)
                h.write(knums+"\t")
        h.write("\n")
        h.close()

def create_tuple_modules(fixed_module_file):
    """
    Generates a tuple from the indication of Modules in which to look for incompleteness, for further use.

    Args:
        fixed_module_file   (str): file name of a ".instruction" file with indication of KEGG Modules of interest

    Returns:
        tuple_modules     (tuple): Python tuple-object including all KEGG Modules of interest
    """
    list_modules = []
    with open(fixed_module_file) as f:
        for line in f.readlines():
            list_modules.append(line.strip())
    tuple_modules = tuple(list_modules)
    return tuple_modules

def create_tuple_modules_1BM(fasta_id, fixed_module_file, oneBM_modules_dir, report_tsv_directory):
    """
    Generates a tuple including Modules missing 1 orthologs block, for further use.

    Args:
        fasta_id                (str): identificative FASTA name for a given MAG/Genome
        fixed_module_file       (str): generic file name of a ".instruction" file with indication of KEGG Modules of interest
        oneBM_modules_dir       (str): output folder of MAG/Genome specific ".instruction" file with KEGG Modules of interest
        report_tsv_directory    (str): testcompleteness_tsv() output folder, to identify 1 block missing modules

    Returns:
        tuple_modules         (tuple): Python tuple-object including all KEGG Modules of interest
    """
    os.chdir(report_tsv_directory)
    list_modules = []
    for file in os.listdir():
        if file.endswith(".tsv") and fasta_id in file:
            with open(file) as f:
                for line in f.readlines():
                    line = line.strip().split("\t")
                    MOD = line[0]
                    COMPLETENESS = line[2]
                    if COMPLETENESS == "1 BLOCK MISSING":
                        list_modules.append(MOD)

    os.chdir(oneBM_modules_dir)
    m = open(fasta_id+"_"+fixed_module_file, "w")
    for module in list_modules:
        m.write(module+"\n")
    m.close()

    tuple_modules = tuple(list_modules)
    return tuple_modules

def write_KOs_from_modules(fasta_id, tuple_modules, report_txt_directory, klists_directory):
    """
    Generates a non-redundant list of KOs to be checked via HMM for Modules of interest,
    either fixed or related to missing annotated genomic content.

    Args:
        fasta_id                    (str): identificative FASTA name for a given MAG/Genome
        tuple_modules             (tuple): output of "create_tuple_modules()" or "create_tuple_modules_1BM()"
        report_txt_directory        (str): testcompleteness() output folder, to identify KOs missing from Modules of interest
        klists_directory            (str): output ".klist" files folder path - in which to save MAG/Genome missing KOs of interest
    """
    os.chdir(report_txt_directory)
    for file in os.listdir():
        if fasta_id in file:
            with open(file) as f:
                klist = []
                v = f.readlines()
                f.seek(0)
                i = 0
                for line in v:
                    i += 1
                    if line.startswith(tuple_modules):
                        start = i
                        j = i
                        for line in v[j:]:
                            j += 1
                            if line.startswith("M0"):
                                end = j
                                break
                        for l in v[start+1:end-2]:
                            l = l.strip().split("\t")[1].split(", ")
                            for KO in l:
                                if not KO in klist:
                                    klist.append(KO)

                os.chdir(klists_directory)
                g = open(file[10:-4]+".klist", "w")
                for KO in klist:
                    g.write(KO+"\n")
                g.close()
                os.chdir(report_txt_directory)

def write_KOs_from_fixed_list(fasta_id, fixed_ko_file, ktests_directory, klists_directory):
    """
    Generates a non-redundant list of KOs to be checked via HMM starting from a fixed list.
    
    Args:
        fasta_id            (str): identificative FASTA name for a given MAG/Genome
        fixed_ko_file       (str): ".instruction" file generated by "setup.py" to be compiled manually with KOs of interest
        ktests_directory    (str): ".ktest" files (KOs file as recovered in pre-annotations) folder path
        klists_directory    (str): output ".klist" files folder path - in which to save MAG/Genome missing KOs of interest
    """
    KO_to_check = []
    os.chdir(dir_base)
    with open(fixed_ko_file) as h:
        for line in h.readlines():
            KO_to_check.append(line.strip())

    os.chdir(ktests_directory)
    for file in os.listdir():
        if fasta_id in file:
            with open(file) as f:
                KO_present = []
                klist = []
                for line in f.readlines():
                    KO = line.strip()
                    KO_present.append(KO)
            for KO in KO_to_check:
                if not KO in KO_present:
                    klist.append(KO)

            os.chdir(klists_directory)
            g = open(file[:-6]+".klist", "w")
            for KO in klist:
                g.write(KO+"\n")
            g.close()
    os.chdir(report_txt_directory)

def taxonomy_filter(taxonomy, dir_base, taxa_file, taxa_dir, update = False):
    """
    Generates a file that includes KEGG Brite species codes (E-level)
    for a given C-level (phylum, most of the times) taxonomy indication.

    Args:
        taxonomy            (str): KEGG Brite taxonomy for MAG/Genome of interest
                                    as indicated in the "genomes.instruction" file.
        dir_base            (str): folder path in which "kemet.py" was executed
        taxa_file           (str): ".keg" output file, that contains each codes of species allowed for subsequent GENES download
        taxa_dir            (str): output folder path
        update   (bool, optional): flag to update KEGG Brite taxonomy - necessary e.g. for the first KEMET execution. Defaults to False.

    Returns:
        taxa_allow         (list): Python-list object including each codes of species allowed for subsequent GENES download
    """
    os.chdir(dir_base)
    taxa_allow=[]
    with open("br08601.keg") as f:
        v = f.readlines()
        f.seek(0)
        i = 0
        for line in v:
            i += 1
            if line.startswith("C") and taxonomy+" " in line:
                i_start = i-1
                break
        for line in v[i_start+1:]:
            i += 1
            if line.startswith("C"):
                i_stop = i-1
                break

        for line in v[i_start:i_stop]:
            if line.startswith("E"):
                taxa_allow.append(line.strip().replace("E        ","").split("  ")[0])

    if update:
        os.chdir(taxa_dir)
        g = open(taxa_file,"w")
        for el in taxa_allow:
            g.write(el+"\n")
        g.close()

    return taxa_allow

def download_ntseq_of_KO(klist_file, dir_base_KO, dir_KO, klists_directory, taxa_dir, taxa_file, base_com_KEGGget):
    """
    Using KEGG API, downloads KEGG flat-files with nt sequences of KOs of interest
    from the allowed species (E-level), following filering of "taxonomy_filter()".

    IF KEGG ACCESS IS AVAILABLE, it is possible to modify "Pool(processes=3)" with processes=N
    and it is possible to download multiple files via API, in full compliance to KEGG license.

    Args:
        klist_file          (str): ".klist" input file name, indicating missing KOs of interest from MAG/Genome
        dir_base_KO         (str): base KEGG KO GENES sequences folder path
        dir_KO              (str): KEGG KO GENES sequences folder path, with taxonomic scope indicated in the command-line input
        klists_directory    (str): ".klist" input files folder path
        taxa_dir            (str): "taxa_file" input file folder path
        taxa_file           (str): ".keg" output file, that contains each codes of species allowed for subsequent GENES download
        base_com_KEGGget    (str): base command of KEGG API "GET" function - which is modified for each entry of GENES
    """
    print(_timeinfo()+"START download nucleotidic sequences")
    os.chdir(taxa_dir)
    taxa_allow = []
    with open(taxa_file) as f:
        for line in f.readlines():
            taxa_allow.append(line.strip())

    os.chdir(klists_directory)
    with open(klist_file) as f:
        os.chdir(dir_base_KO)
        if not os.path.exists(dir_KO):
            os.mkdir(dir_KO)
        os.chdir(dir_KO)
        for line in f.readlines():
            line = line.strip()
            flatfile = str(line)+".keg"
            cmd_dir = "mkdir "+line
            os.chdir(dir_KO)
            if not os.path.exists(dir_KO+line):
                os.system(cmd_dir)
            else:
                continue
            os.chdir(line)
            os.system(base_com_KEGGget+line+" > "+flatfile)

            genes = parsekoflat(flatfile)
            os.system("rm "+flatfile)

            if __name__ == '__main__':
                # requests to KEGG API without a granted access are limited (check KEGG LICENCE)
                # POSSIBILITY: modify next line "(processes= n)" if access to KEGG is available
                with Pool(processes=3) as p:
                    p.map(getntseq, genes)
                p.close()
    print(_timeinfo()+"COMPLETE download nucleotidic sequences")

def parsekoflat(file):
    """
    Parses KO flatfiles obtained from KEGG API,
    in order to generate the filtered list of sequences for a bulk download.
    (Called within the "download_ntseq_of_KO()" function).

    Args:
        file        (str): KEGG API KO flat-file file name

    Returns:
        genes      (list): Python-list object with genes connected to the KO,
                            for each appropriate species within the specified BRITE taxonomy
    """
    genes = []
    with open(file) as f:
        v = f.readlines()
        n = 0
        for line in v:
            if not line.startswith("GENES"):
                n += 1
            else:
                break
        f.seek(0)
        for line in v[n:]:
            if line.startswith("REFERENCE") or line.startswith("///"):
                break
            line = line.replace("GENES       ","").strip()
            m = line.index(":")
            species = line[:m+1].casefold()
            line_s = line.split()
            for g_name in line_s[1:]:
                gene = (species+g_name).casefold()
                if "(" in gene:
                    p = gene.index("(")
                    gene = gene[:p]
                if "draft" in gene:
                    continue
                genes.append(gene)
    return genes

def getntseq(gene):
    """
    Downloads nt sequence of a given gene, from the list of KO-related genes list.
    (Called within the "download_ntseq_of_KO()" function).

    Args:
        gene        (str): element of "genes" input list, generated via "parsekoflat()"

    Returns:
        True       (bool): only necessary for script continuation
    """
    gene_name = gene
    stop = gene_name.find(":")
    gene_taxa = gene_name[:stop]

    if gene_taxa in taxa_allow:
        cmd_get_ntseq = base_com_KEGGget+gene_name+"/ntseq"
        os.system(cmd_get_ntseq+" > "+gene_name+".fna")
        return 1

def filter_and_allign(taxa_dir, taxa_file, fasta_id, klist_file, klists_directory, msa_dir, dir_KO):
    """
    Generates a nucleotidic multifasta with sequences from the given taxonomy range.
    The output does NOT contain redundant sequences.
    Keeps results in a folder organized by the FASTA id of MAG/Genome.

    Args:
        taxa_dir            (str): "taxa_file" input file folder path
        taxa_file           (str): ".keg" file, contains each codes of species allowed for subsequent GENES download
        fasta_id            (str): identificative FASTA name for a given MAG/Genome
        klist_file          (str): ".klist" file name, indicating missing KOs of interest from MAG/Genome
        klists_directory    (str): ".klist" input files folder path
        msa_dir             (str): ".fna" nt multifasta (single representative sequences) output folder path
        dir_KO              (str): KEGG KO GENES sequences folder path, with taxonomic scope indicated in the command-line input
    """
    print(_timeinfo()+"START sequences filtering and allignment")
    ### filter for taxa of interest
    os.chdir(taxa_dir)
    taxa_allow = []
    with open(taxa_file) as f:
        for line in f.readlines():
            taxa_allow.append(line.strip())

    os.chdir(msa_dir)
    if not fasta_id in os.listdir():
        os.mkdir(fasta_id)

    os.chdir(klists_directory)
    KO_to_align = []
    with open(klist_file) as f:
        for line in f.readlines():
            KO = line.strip()
            if not KO in KO_to_align:
                KO_to_align.append(KO)

    os.chdir(dir_KO)
    for K in os.listdir():
        if not K in KO_to_align:
            continue
    # dictionary of non-redundant nt sequences (100% identity)
    # in order not to overvalue species with different strains in KEGG taxonomy
    # but only focusing on SEQUENCE DIVERSITY
        os.chdir("./"+K)
        sequniq = {} # {sequence : tax_code_of_identical_seqs}
        for nt_file in os.listdir():
            code = nt_file.split(":")[0]
            if not code in taxa_allow:
                continue
    ### exclude redundant nt sequences
            with open(nt_file) as f:
                seq = f.readlines()[1:]
                seq1 = "".join(seq).replace("\n","")

                if not seq1 in sequniq.keys():
                    vett = [nt_file]
                    sequniq.update({seq1:vett})
                else:
                    vett = sequniq[seq1]
                    vett.append(nt_file)
                    sequniq.update({seq1:vett})
    ### Write a multiple sequence fasta
        os.chdir(msa_dir+fasta_id)
        if not K in os.listdir():
            os.mkdir(K)
        os.chdir(msa_dir+fasta_id+"/"+K)
        f = open("MSA_"+K+".fna", "a")
        for key, value in sequniq.items():
            f.write(">"+str(value[0][:-4])+"\n")
            f.write(key+"\n")
        f.close()
        os.chdir(dir_KO)
    print(_timeinfo()+"COMPLETE Filter and allign")

def MSA_and_HMM(msa_dir_comm, base_com_mafft, base_com_hmmbuild, log=False):
    """
    Runs MAFFT alignment and then build nt profile HMM from it.
    Keeps results in a folder organized by the FASTA-header of MAG/Genome.

    Args:
        msa_dir_comm                    (str): nt multi-fasta folder path, as modified for the MAG/Genome of interest
        base_com_mafft                  (str): base command for MAFFT execution - modified for each entry of GENES (KO)
        base_com_hmmbuild               (str): base command for "hmmbuild" execution - modified for each entry of GENES (KO)
        log                  (bool, optional): keep execution times in a log file (if specified in command-line args). Defaults to False.
    """
    print(_timeinfo()+"START MSA and HMMs creation")
    os.chdir(msa_dir_comm)
    for K in os.listdir():
        os.chdir(K)
        ch_com_mafft = base_com_mafft.replace("K_NUMBER", K)
        ch_com_hmmbuild = base_com_hmmbuild.replace("K_NUMBER", K)
        os.system(ch_com_mafft)
        if log==True:
            logging.info('COMPLETE MAFFT execution')
        os.system(ch_com_hmmbuild)
        if log==True:
            logging.info('COMPLETE hmmbuild execution')
        os.chdir(msa_dir_comm)
    print(_timeinfo()+"COMPLETE MSA and HMM creation")

def nhmmer_for_genome(fasta_genome, msa_dir_comm, base_com_nhmmer):
    """
    Runs a nHMMER search for newly generated HMM profiles against a given MAG/Genome.
    Keeps results in a folder organized by the FASTA id of MAG/Genome.

    Args:
        fasta_genome        (str): identificative FASTA name (including path) for a given MAG/Genome
        msa_dir_comm        (str): nt multi-fasta folder path, as modified for the MAG/Genome of interest
        base_com_nhmmer     (str): base command for "nhmmer" execution - modified for each entry of GENES (KO)
    """
    print(_timeinfo()+"START nhmmer search")
    os.chdir(msa_dir_comm)
    for K in os.listdir():
        os.chdir(K)
        ch_com_nhmmer = base_com_nhmmer.replace("K_NUMBER", K).replace("PATHFILE", fasta_genome)
        os.system(ch_com_nhmmer)
        os.chdir(msa_dir_comm)
    print(_timeinfo()+"COMPLETE nhmmer")

def move_HMM_and_clean(hmm_dir_comm, msa_dir_comm):
    """
    Orders HMMs, moving them from MSA folder into a dedicated HMM folder.
    Removes multi-fasta sequences (".fna") necessary for MAFFT proper alignment.

    Args:
        hmm_dir_comm    (str): HMM folder path, as modified for the MAG/Genome of interest
        msa_dir_comm    (str): nt multi-fasta folder path, as modified for the MAG/Genome of interest
    """
    os.chdir(hmm_dir)
    if not os.path.exists(hmm_dir_comm):
        os.mkdir(hmm_dir_comm)
    os.chdir(msa_dir_comm)
    K_numbers = os.listdir()
    for K in K_numbers:
        os.chdir(hmm_dir_comm)
        if not os.path.exists(hmm_dir_comm+K):
            os.mkdir(K)
        os.chdir(msa_dir_comm+K)
        for file in os.listdir():
            if file.endswith(".hmm") or file.endswith(".hits"):
                hmm_file = file
                os.replace(msa_dir_comm+K+"/"+hmm_file, hmm_dir_comm+K+"/"+hmm_file)
            if file.endswith(".fna"):
                os.remove(file)
    print(_timeinfo()+"COMPLETE move HMM and clean")

def movebackHMM(hmm_dir_comm, msa_dir_comm):
    """
    Moves ".hmm" files back to multi-alignment folder,
    e.g. in order to try different HMM-hits threshold.
    WHEN TO USE THIS:
    to unroll from a complete command-line run - called with the arg: "--retry_nhmmer"

    Args:
        hmm_dir_comm    (str): HMM folder path, as modified for the MAG/Genome of interest
        msa_dir_comm    (str): nt multi-fasta folder path, as modified for the MAG/Genome of interest
    """
    os.chdir(hmm_dir_comm)
    K_numbers = os.listdir()
    for K in K_numbers:
        os.chdir(hmm_dir_comm+K)
        for file in os.listdir():
            if file.endswith(".hmm"):
                hmm_file = file
                os.replace(hmm_dir_comm+K+"/"+hmm_file, msa_dir_comm+K+"/"+hmm_file)
    print(_timeinfo()+"COMPLETE move back HMM")

def nhmmer_significant_hits_corr(fasta_id, hmm_dir_comm, threshold= 100, corr_threshold= _def_thr, evalue_threshold = float(1e-30)):
    """
    Stores the nhmmer-derived HMM hits in a dedicated ".txt" file.
    Detailed informations stored include the scoring, genomic context and
    HMM positions of the hits found with nhmmer.

    Args:
        fasta_id                      (str): identificative FASTA name for a given MAG/Genome
        hmm_dir_comm                  (str): HMM folder path, as modified for the MAG/Genome of interest
        threshold           (int, optional): Minimum nhmmer score for a hit to be considered proper. Defaults to 100.
        corr_threshold    (float, optional): Minimum nhmmer score, corrected by profile HMM lenght
                                             for a hit to be considered proper. Defaults to _def_thr.
        evalue_threshold  (float, optional): Minimum nhmmer e-value for a hit to be considered proper
                                             (not in use in current form). Defaults to float(1e-30).

    Returns:
        sig_hits                     (dict): Python-dictionary object including HMM results
                                             (not strictly necessary for script continuation).
    """
    os.chdir(hmm_dir_comm)
    sig_hits = {}

    g = open(fasta_id+"_HMM_hits.txt","a")
    g.write("K\tcorr_score,evalue\tfragment\tstrand\tl_bound\tr_bound\tp_lenght\thmmfrom\thmmto\n")
    g.close()

    for directory in sorted(os.listdir()):
        if directory.startswith("K"):
            K = directory
            os.chdir("./"+K)
            if not K+".hits" in os.listdir():
                os.chdir(hmm_dir_comm)
            else:
                evalue = 100
                score = 0
                corr_score = 0
                with open(K+".hits") as f:
                    v = f.readlines()[:-9]          # NOT INCLUDING lines w/ the nhmmer command as they rise an error
                    add_ = int(v[1].index("-")-1)   # CORRECT for long contig names
                    for line in v:
                        if K in line[32+add_:38+add_]:
                            evalue = line[127+add_:136+add_].strip()
                            score = line[137+add_:143+add_].strip()
                            fragment = line[0:20+add_].strip()
                            strand = line[120+add_:126+add_].strip()
                            left_bound = int(line[80+add_:87+add_].strip())
                            right_bound = int(line[88+add_:95+add_].strip())
                            bounds = [str(left_bound-1),str(right_bound-1)]
                            hmmfrom = int(line[65+add_:71+add_].strip())
                            hmmto = int(line[72+add_:79+add_].strip())
                            profile_lenght = int(hmmto - hmmfrom)

                            corr_score = round((float(score)/profile_lenght), 4)
                            break

                    if corr_score > corr_threshold and float(score) > threshold: # could also be unified with evalue_threshold
                        sig_hits.update({K:[fragment,strand,bounds]})
                        os.chdir(hmm_dir_comm)
                        g = open(fasta_id+"_HMM_hits.txt","a")
                        g.write(K+"\t"+str(corr_score)+","+str(evalue)+"\t"+fragment+"\t"+strand+"\t"+str(left_bound-1)+"\t"+str(right_bound-1)+"\t"+str(profile_lenght)+"\t"+str(hmmfrom-1)+"\t"+str(hmmto-1)+"\n")
                        g.close()
                os.chdir(hmm_dir_comm)

    os.rename(hmm_dir_comm+fasta_id+"_HMM_hits.txt",hmm_hits_dir+fasta_id+"_HMM_hits.txt")
    print(_timeinfo()+"COMPLETE nhmmer significant hits")
    return sig_hits

def HMM_hits_sequences(hmm_hits_dir, dir_genomes):
    """
    Loads MAGs/Genomes HMM hits info, previously generated via "nhmmer_significant_hits_corr()".
    Checks hits in genome content (contigs-scaffolds-chromosomes), both in FOR & REV directions.
    Stores relevant hits in a dictionary-object.

    Args:
        hmm_hits_dir        (str): HMM hits folder path
        dir_genomes         (str): command-line defined genomes folder path

    Returns:
        HMM_hits_dict      (dict): Python-dictionary object:
                                    keys: ">MAG/Genome+KO"
                                    values: "nt sequence"
    """
    #### load all infos of different genomes' HMM-hits
    MAG_Khit_dict = {}
    HMM_hits_dict = {}

    os.chdir(hmm_hits_dir)
    for el in sorted(os.listdir()):
        if not el.endswith(".txt"):
            continue
        MAG = el[:-13] # slice for "_HMM_hits.txt"
        v_Khits = []
        with open(el) as f:
            for line in f.readlines()[1:]:
                line = line.strip().split("\t")
                v_Khit_info = []
                K = line[0]
                fragment = line[2]
                strand = line[3]
                l_bound = line[4]
                r_bound = line[5]

                v_Khit_info.append(K)
                v_Khit_info.append(fragment)
                v_Khit_info.append(strand)
                v_Khit_info.append(l_bound)
                v_Khit_info.append(r_bound)
                v_Khits.append(v_Khit_info)
        MAG_Khit_dict.update({MAG:v_Khits})

    #### once all infos of the hits are loaded, check the contigs//genomes for hits sequences
    for genome, val in MAG_Khit_dict.items():
        for hits in MAG_Khit_dict[genome]:
            strand_plus = False
            K = hits[0]
            fragment = hits[1]
            strand = hits[2]
            l_bound = int(hits[3])
            r_bound = int(hits[4])

            os.chdir(dir_genomes)
            for genome_file in os.listdir():
                if genome_file.startswith(genome+"."):
                    with open(genome_file) as f:
                        v_lines = f.readlines()
                        f.seek(0)
            i = -1
            for line in v_lines:
                line = line.strip()
                i += 1
                if ">" and fragment in line:
                    i_frag_start = i
                    j = i
                    for line in v_lines[i_frag_start+1:]:
                        j += 1
                        if ">" in line or j >= (len(v_lines)-1): # new contigs or file-end
                            i_frag_stop = j
                            break
                    break
    #### load contigs' "+" strand sequence
            try:
                fragment_w_hit = ""
                for line in v_lines[i_frag_start+1:i_frag_stop]:
                    line = line.strip()
                    for char in line:
                        fragment_w_hit += char
    #### different treatment for hits in "+" & "-" strands
    #### do not include hits with "N", resulting e.g. from scaffold placeholders
                if strand == "+":
                    strand_plus = True
                if strand_plus == True:
                    SEQUENCE = fragment_w_hit[l_bound:r_bound]
                    if not "N" in SEQUENCE and not "n" in SEQUENCE:
                        HMM_hits_dict.update({">"+genome+"_"+K:SEQUENCE})
                elif strand_plus == False:
                    SEQUENCE_pre = fragment_w_hit[r_bound:l_bound]
                    pre = "ACTG"
                    post = "TGAC"
                    compl = SEQUENCE_pre.maketrans(pre, post)
                    seq_compl = SEQUENCE_pre.translate(compl)
                    SEQUENCE = seq_compl[::-1]
                    if not "N" in SEQUENCE and not "n" in SEQUENCE:
                        HMM_hits_dict.update({">"+genome+"_"+K:SEQUENCE})
            except:
                pass
    return HMM_hits_dict

def HMM_hits_translated_sequences(HMM_hits_dict):
    """
    Translates HMM hits in the 3 frames of the adequate sequence strand.
    Stores them in a dictionary-object.

    Args:
        HMM_hits_dict               (dict): Python-dictionary object output of "HMM_hits_sequences()"

    Returns:
        HMM_hits_TRANSLATED_dict    (dict): Python-dictionary object:
                                            keys: ">MAG/Genome+KO+frame_indication"
                                            values: "aa sequence"
    """
    #POSSIBILITY: NOT ONLY t11 table, but also others
    t11 = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
           "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
           "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
           "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
           "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
           "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
           "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
           "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
           "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
           "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
           "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
           "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
           "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
           "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
           "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
           "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

    HMM_hits_TRANSLATED_dict = {}
    for fasta_id, seq in HMM_hits_dict.items():
        ntseq1 = seq.upper()
        ntseq2 = seq[1:].upper()
        ntseq3 = seq[2:].upper()
        fasta1 = fasta_id+"__f1"
        fasta2 = fasta_id+"__f2"
        fasta3 = fasta_id+"__f3"
        aaseq1 = ""
        aaseq2 = ""
        aaseq3 = ""

        codons1 = [ntseq1[i:i+3] for i in range(0, len(ntseq1), 3)]
        for el in codons1:
            if len(el) < 3:
                codons1.remove(el)
        for codon in codons1:
            tr_codon = t11[codon]
            aaseq1 += tr_codon
        HMM_hits_TRANSLATED_dict.update({fasta1:aaseq1})

        codons2 = [ntseq2[i:i+3] for i in range(0, len(ntseq2), 3)]
        for el in codons2:
            if len(el) < 3:
                codons2.remove(el)
        for codon in codons2:
            tr_codon = t11[codon]
            aaseq2 += tr_codon
        HMM_hits_TRANSLATED_dict.update({fasta2:aaseq2})

        codons3 = [ntseq3[i:i+3] for i in range(0, len(ntseq3), 3)]
        for el in codons3:
            if len(el) < 3:
                codons3.remove(el)
        for codon in codons3:
            tr_codon = t11[codon]
            aaseq3 += tr_codon
        HMM_hits_TRANSLATED_dict.update({fasta3:aaseq3})

    return HMM_hits_TRANSLATED_dict

def HMM_hits_longest_translated_sequences(HMM_hits_dict, HMM_hits_TRANSLATED_dict):
    """
    Compares the HMM translated hits and keep the longest hit without a stop codon.
    Stores them in a dictionary-object.

    Args:
        HMM_hits_dict                   (dict): Python-dictionary object output of "HMM_hits_sequences()"
        HMM_hits_TRANSLATED_dict        (dict): Python-dictionary object output of "HMM_hits_translated_sequences()"

    Returns:
        HMM_hits_TRANSLATED_MAXLEN_dict (dict): Python-dictionary object:
                                                keys: ">MAG/Genome+KO+frame_indication"
                                                values: "aa sequence" - only the longest without stop codons (*)
    """
    max_len_dict = {}
    for fasta_nf in HMM_hits_dict.keys():
        max_len_dict.update({fasta_nf:[]})

    HMM_hits_TRANSLATED_MAXLEN_dict = {}

    for fasta_id, seq in HMM_hits_TRANSLATED_dict.items():
        fasta_nf = fasta_id.split("__")[0]
        seq = seq.split("*")                   # divided by stop codons
        seq_max = max(seq, key = len)          # longer for single reading frame
        max_len_dict[fasta_nf].append(seq_max) # add the longest to list

    for fasta_nf in max_len_dict.keys():
        seq_max_allframes = max(max_len_dict[fasta_nf], key = len)
        frame = "__f"+str(max_len_dict[fasta_nf].index(seq_max_allframes)+1)
        HMM_hits_TRANSLATED_MAXLEN_dict.update({fasta_nf+frame:seq_max_allframes})

    return HMM_hits_TRANSLATED_MAXLEN_dict

def recap_hits_corr(fasta_id, hmm_hits_dir, HMM_hits_dict, HMM_hits_longestTRANSLATED_dict, run_start):
    """
    Stores informations from the HMM-hits dictionaries, alongside identified sequences.
    It gives the final tabular output of the entire HMM procedure.

    Args:
        fasta_id                            (str): identificative FASTA name for a given MAG/Genome
        hmm_hits_dir                        (str): HMM hits folder path
        HMM_hits_dict                      (dict): Python-dictionary object output of "HMM_hits_sequences()"
        HMM_hits_longestTRANSLATED_dict    (dict): Python-dictionary object output of "HMM_hits_longest_translated_sequences()"
        run_start                           (str): indication of the date in which the "kemet.py" process started
    """
    os.chdir(hmm_hits_dir)
    current_run = "KEMET_run_"+run_start
    if not current_run in os.listdir():
        os.mkdir(current_run)
    os.chdir(current_run)
    f = open("file_recap_"+run_start+".tsv", "a")
    header = "MAG\tKO\tcorr_score,evalue\tfragment\tstrand\tl_bound\tr_bound\tp_lenght\thmmfrom\thmmto\tframe\tseq\txseq\n"
    if f.tell() == 0: 
        f.write(header)

    os.chdir(hmm_hits_dir)
    for file in sorted(os.listdir()):
        if file.endswith(".txt") and fasta_id in file:
            with open(file) as g:
                for line in g.readlines()[1:]:
                    KO = line.strip().split("\t")[0]
                    for hit in HMM_hits_dict.keys():
                        if hit.startswith(">"+file[:-13]+"_"+KO): ###
                            seq = HMM_hits_dict[hit].upper()
                            break
                    for hit in HMM_hits_longestTRANSLATED_dict.keys():
                        if hit.startswith(">"+file[:-13]+"_"+KO): ###
                            frame = hit[-1]
                            xseq = HMM_hits_longestTRANSLATED_dict[hit]
                            break
                    os.chdir(current_run)
                    f.write(file[:-13]+"\t"+line.strip()+"\t"+frame+"\t"+seq+"\t"+xseq+"\n")
                    os.chdir(hmm_hits_dir)
    f.close()

def build_de_novo_GSMM(FASTA, fasta_genome, de_novo_model_directory, current_run, log=False):
    """
    Generates de-novo GSMM from a contextual Prodigal gene calling,
    adding HMM-hits derived translated hits to the predicted proteome (".faa").

    Args:
        FASTA                           (str): identificative FASTA name for a given MAG/Genome - with extension
        fasta_genome                    (str): identificative FASTA name (including path) for a given MAG/Genome
        de_novo_model_directory         (str): output files folder path
        current_run                     (str): indication of the date in which the "kemet.py" process started
        log                  (bool, optional): keep execution times in a log file (if specified in command-line args). Defaults to False.
    """
    fasta_id = FASTA.replace(".fasta","").replace(".fna","").replace(".fa","")
    HMM_HITS = {}
    os.chdir(hmm_hits_dir+current_run)
    for file in os.listdir():
        if not file.startswith("file_recap_") and file.endswith(".tsv"): # open the hmm recap
            continue
        with open(file) as f:
            for line in f.readlines()[1:]:
                line = line.strip().split("\t")
                MAG = line[0]
                KO = line[1]
                xseq = line[12]
                if not MAG == fasta_id:
                    continue
                HMM_HITS.update({">"+KO:xseq})

    # launch prodigal gene calling - default:QUIET
    os.chdir(de_novo_model_directory)
    if not "proteins" in os.listdir():
        os.mkdir("proteins")
    prodigal_cmd = "prodigal -i fasta_genome -a ./proteins/FASTA.faa -q -m > /dev/null"
    prodigal_cmd_mod = prodigal_cmd.replace("fasta_genome",fasta_genome).replace("FASTA",fasta_id)
    os.system(prodigal_cmd_mod)
    print(_timeinfo()+"COMPLETE Prodigal command for "+FASTA)
    if log == True:
        logging.info('COMPLETE Prodigal command for '+FASTA)

    # ADD HMM-hits
    os.chdir(de_novo_model_directory+"proteins")
    f = open(fasta_id+".faa", "a")
    for ko, xseq in HMM_HITS.items():
        f.write(ko+"\n")
        f.write(xseq+"\n")
    f.close()
    os.chdir(de_novo_model_directory)

    # GENERATE CARVEME MODEL - default:QUIET
    if not "dmnd_intermediates" in os.listdir():
        os.mkdir("dmnd_intermediates")
    carveme_cmd = "carve ./proteins/FASTA.faa --fbc2 -u universe -o FASTA.xml 2> /dev/null"
    carveme_cmd_mod = carveme_cmd.replace("FASTA",fasta_id).replace("universe",metabolic_universe)
    print(_timeinfo()+"START CarveMe command for "+FASTA)
    if log == True:
        logging.info('START CarveMe command for '+FASTA)
    os.system(carveme_cmd_mod)
    print(_timeinfo()+"COMPLETE CarveMe command for "+FASTA)
    if log == True:
        logging.info('COMPLETE CarveMe command for '+FASTA)

    os.chdir(de_novo_model_directory+"/proteins/")
    move_dmnd_cmd = "mv "+fasta_id+".tsv ../dmnd_intermediates"
    os.system(move_dmnd_cmd)

def list_all_modules(Modules_directory):
    """
    Helper function to list all KEGG Modules available
    for KEMET Genome-scale models procedures.

    Args:
        Modules_directory   (str): KEGG Modules flat-files folder path

    Returns:
        list_all_mod       (list): Python-list object including all KEGG Modules available for further use
    """
    os.chdir(Modules_directory)
    list_all_mod = []
    all_mod = sorted(os.listdir())
    for module in all_mod:
        if module.endswith(".txt"):
            list_all_mod.append(module)
    return list_all_mod

def searchKeggShort(list_all_mod, Modules_directory, knumber):
    """
    Returns a list-object of Module ids
    format: [Mxxxxx], in which a given input KO is found.

    Args:
        list_all_mod           (list): output of "list_all_modules()" function
        Modules_directory       (str): KEGG Modules flat-files folder path
        knumber                 (str): input KO to be searched in all KEGG Modules definitions

    Returns:
        hits                   (list): Python-list object including all KEGG Modules in which input KO is present
    """
    os.chdir(Modules_directory)
    hits = []
    for element in list_all_mod:
        name = element
        with open(name) as f:
            n = 0
            manydefinitions = []
            for line in f.readlines()[2:]:
                n += 1
                if not line.startswith("ORTHOLOGY") and n == 1:
                    definition = line.replace("-- ", "").strip().split("  ", 1)[1]
                    continue
                if line.startswith("ORTHOLOGY"):
                    if len(manydefinitions) != 0:
                        definition = manydefinitions
                    break
                if not line.startswith("ORTHOLOGY") and n > 1:
                    manydefinitions.append(definition)
                    definition = line.strip()
                    manydefinitions.append(definition)
                    continue

        if len(manydefinitions) > 0:
            for singledef in manydefinitions:
                if knumber in singledef:
                    hits.append(name[:6])
                    break
                else:
                    continue
            else:
                continue
        elif len(manydefinitions) == 0:
            if knumber in definition:
                hits.append(name[:6])
                continue
            else:
                continue
    if len(hits) == 0:
        return knumber+" no hits!" # RAISE error?
    else:
        return hits

def knum4reac_mod(file):
    """
    Returns a dictionary-object that connects KOs to the reactions (R),
    as identified from KEGG Module flat-files as keys and reactions (R) as values, as per Module file.
    (Called in the "total_R_from_KOhits()" function).
    
    Args:
        file                        (str): input KEGG Module flat-file path

    Returns:
        dictionary_knumber_reac    (dict): Python-dictionary object:
                                            keys: "KO" entry
                                            values: "R" entry
    """
    with open(file) as f:
        dictionary_knumber_reac = {}
        l_count = 0
        for line in f.readlines():
            if not line.startswith("ORTHOLOGY"):
                l_count += 1
                continue
            if line.startswith("ORTHOLOGY"):
                break

        f.seek(0)
        for line in f.readlines()[l_count:]:
            if line.startswith("CLASS"):
                break
            if line.startswith("ORTHOLOGY"):
                line = line.replace("ORTHOLOGY", "")
            if not "[RN:" in line:
                continue
            knumber_a = line.strip().split("  ")[0]
            knumber = re.split("[+-,]", knumber_a)
            reaction = line.strip().split("[RN:")[1].replace("]", "").replace(" ",",")
            for knum_element in knumber:
                if not knum_element in dictionary_knumber_reac.keys():
                    dictionary_knumber_reac.update({knum_element:reaction})
                else:
                    prior_reaction = str(dictionary_knumber_reac[knum_element])
                    new_reaction = prior_reaction+","+reaction
                    dictionary_knumber_reac.update({knum_element:new_reaction})
    f.close()

    return dictionary_knumber_reac

def KEGG_BiGG_SEED_RN_dict(reactions_DB, DB_directory, ontology = "BiGG"):
    """
    Connects KEGG RN to BiGG IDs or KEGG RN to ModelSEED IDs
    starting from a tabular data file which connects ModelSEED reactions
    to other GSMM reaction onthologies.

    Args:
        reactions_DB                (str): input ".tsv" reaction file path
        DB_directory                (str): input file folder path
        ontology          (str, optional): Genome-scale models ontology indication. Defaults to "BiGG".

    Returns:
        dict_kegg_x_dict           (dict): Python-dictionary object with manyXmany relations of KEGG Rs and BiGG IDs
        dict_kegg_x_seed           (dict): (Unsupported in current version)
    """
    dict_kegg_x_dict = {}
    dict_kegg_x_seed = {}
    os.chdir(DB_directory)
    with open(reactions_DB) as g:
        for line in g.readlines()[1:]:
            x = line.split("\t")
            #variable names
            rxn = x[0]
            bigg_ids = x[2:24]
            kegg_ids = x[25:34]

            #filter for empty lists
            non_empty_bigg_ids = list(filter(None, bigg_ids))
            non_empty_kegg_ids = list(filter(None, kegg_ids))

            #unique lists
            unique_kegg_ids = []
            unique_bigg_ids = []

            #dict added for every reaction/line
            dict_rxn_x_bigg = {}
            #only non-redundant KEGG RN entry x reaction/line
            for a in non_empty_kegg_ids:
                if not a in unique_kegg_ids:
                    unique_kegg_ids.append(a)

            #only non-redundant BiGG entry x reaction/line
            for b in non_empty_bigg_ids:
                if not b in unique_bigg_ids:
                    unique_bigg_ids.append(b)

        # {KEGG : [{rxns:[BiGGs]}]}
            #for each KEGG RN, create-update a non-redundant dict-entry
            for Kegg in unique_kegg_ids:
                if not Kegg in dict_kegg_x_dict.keys():
                    vett_dict_rxn_bigg=[]
                    dict_rxn_x_bigg.update({rxn:unique_bigg_ids})
                    vett_dict_rxn_bigg.append(dict_rxn_x_bigg)
                    dict_kegg_x_dict.update({Kegg:vett_dict_rxn_bigg})
                else: #Kegg in dizio_kegg_x_diz.keys():
                    vett_dict_rxn_bigg=dict_kegg_x_dict[Kegg]
                    dict_rxn_x_bigg.update({rxn:unique_bigg_ids})
                    vett_dict_rxn_bigg.append(dict_rxn_x_bigg)
                    dict_kegg_x_dict.update({Kegg:vett_dict_rxn_bigg})

        # {KEGG : [rxns]}
            #for each KEGG RN, create-update a non-redundant dict-entry
            for Kegg in unique_kegg_ids:
                if not Kegg in dict_kegg_x_seed.keys():
                    vett_dict_rxn_seed=[]
                    vett_dict_rxn_seed.append(rxn)
                    dict_kegg_x_seed.update({Kegg:vett_dict_rxn_seed})
                else: #Kegg in dict_kegg_x_seed.keys():
                    vett_dict_rxn_seed=dict_kegg_x_seed[Kegg]
                    vett_dict_rxn_seed.append(rxn)
                    dict_kegg_x_seed.update({Kegg:vett_dict_rxn_seed})

    if ontology == "BiGG":
        return dict_kegg_x_dict
    elif ontology == "SEED":
        return dict_kegg_x_seed

def keggR_in_DB(list_of_unique_Rnumbers, DB_Kegg_reactions):
    """
    Helper function to search for unique KEGG REACTIONS (R) in
    the reaction database file provided in the KEMET GitHub.

    Args:
        list_of_unique_Rnumbers     (list): Python-list output of "total_R_from_KOhits()" function
        DB_Kegg_reactions           (dict): Python-dictionary output of "KEGG_BiGG_SEED_RN_dict()" function

    Returns:
        v_KeggR_unique_in_DB        (list): Python-list object including all different and unique KEGG R
                                            as recovered from the reaction_DB file
    """
    v_KeggR_unique = []
    v_KeggR_unique_in_DB = []
    for Rnumber in list_of_unique_Rnumbers:
        v_KeggR_unique.append(Rnumber)
        if Rnumber in DB_Kegg_reactions.keys():
            v_KeggR_unique_in_DB.append(Rnumber)

    return v_KeggR_unique_in_DB

def bigg_nonredundant(v_KeggR_unique_in_DB, DB_Kegg_reactions):
    """
    Helper function to search for unique BiGG REACTIONS in
    the reaction database file provided in the KEMET GitHub.

    Args:
        v_KeggR_unique_in_DB    (list): Python-list output of "keggR_in_DB()" function
        DB_Kegg_reactions       (dict): Python-dictionary output of "KEGG_BiGG_SEED_RN_dict()" function

    Returns:
        v_bigg_nonredundant     (list): Python-list object including unique BiGG reaction,
                                        as linked to KEGG R via BiGG database
    """
    v_bigg_nonredundant = []
    for Rnumber in v_KeggR_unique_in_DB:
        for subdict in DB_Kegg_reactions[Rnumber]:
            for v_bigg in subdict.values():
                for single_bigg in v_bigg:
                    if single_bigg not in v_bigg_nonredundant:
                        v_bigg_nonredundant.append(single_bigg)

    return v_bigg_nonredundant

def modelseed_nonredundant(v_KeggR_unique_in_DB, DB_Kegg_reactions):
    """
    Helper function to search for unique ModelSEED REACTIONS in
    the reaction database file provided in the KEMET GitHub.

    Args:
        v_KeggR_unique_in_DB        (list): Python-list output of "keggR_in_DB()" function
        DB_Kegg_reactions           (dict): Python-dictionary output of "KEGG_BiGG_SEED_RN_dict()" function

    Returns:
        v_modelseed_nonredundant    (list): Python-list object including unique BiGG reaction,
                                            as linked to ModelSEED reactions via BiGG database
    """
    v_modelseed_nonredundant = []
    for Rnumber in v_KeggR_unique_in_DB:
        for single_rxn in DB_Kegg_reactions[Rnumber]:
            if not single_rxn in v_modelseed_nonredundant:
                v_modelseed_nonredundant.append(single_rxn)

    return v_modelseed_nonredundant

def bigg_gapfill_absent_in_model(bigg_nonredundant_file, model):
    """
    Helper function to identify putative gap-fill reactions
    missing in the input existing ".xml" model.
    (Called in the "reframed_reaction_addition()" function).

    Args:
        bigg_nonredundant_file              (str): log ".txt" file generated via "log_bigg_nr()"
        model                           (CBModel): Python-object of an existing ".xml" model, as loaded via Reframed

    Returns:
        bigg_gapfill_absent_in_model       (list): Python-list object including unique BiGG reaction,
                                                    not present in the input ".xml" genome-scale model
    """
    biggreactions_names = []
    putative_gapfill = []
    bigg_gapfill_absent_in_model = []

    for reaction in model.reactions.keys():
        if reaction == 'Growth' or reaction == 'R_ATPM':
            continue
        nome_bigg = reaction[2:]
        biggreactions_names.append(nome_bigg)

    with open(bigg_nonredundant_file) as f:
        for line in f.readlines():
            line_s = line.strip()
            putative_gapfill.append(line_s)

    for bigg1 in putative_gapfill:
        if not bigg1 in biggreactions_names:
            bigg_gapfill_absent_in_model.append(bigg1)

    return bigg_gapfill_absent_in_model

def curl_bigg_reaction(single_reac):
    """
    Helper function to download a file that contains BiGG reaction,
    stores it in an intermediate file (".ift").
    (Called in the "reframed_reaction_addition()" function).

    Args:
        single_reac     (str): name of a single reaction, taken from the list-output of "bigg_gapfill_absent_in_model()"

    Returns:
        file_ift        (str): name of intermediate flat file (".ift") including a reaction string
    """
    cmd_line = "curl --silent 'http://bigg.ucsd.edu/api/v2/universal/reactions/"
    el = single_reac
    command = cmd_line+el+"'"
    file_ift = str(el)+".ift"
    os.system(command+" > "+file_ift)

    return file_ift

def api_file_reorder(file_ift, verbose = False):
    """
    Helper function to clean and reorganize info from
    intermediate files, containing reaction strings.
    (Called in the "reframed_reaction_addition()" function).

    Args:
        file_ift            (str): name of intermediate flat file (".ift") including a reaction string
        verbose  (bool, optional): print more info regarding process status. Defaults to False.
    """
    with open(file_ift) as f:
        for line in f.readlines():
            text = "".join(line)
        #text cleaning
            text2 = text.replace(",", ",\n").replace(":",":\t").replace("&#8652;","<->")
            file_txt = str(file_ift).replace(".ift", ".txt")
            #files_txt.append(file_txt)
            g = open(file_txt, "w")
            g.write(text2)
            g.close()
            os.system("rm "+file_ift)

            if verbose:
                print("Info from file "+file_ift+" reordered")
    return

def api_file_reorder2(file_ift, verbose = False):
    """
    Helper function #2 to clean and reorganize info from
    intermediate files, containing reaction strings.

    Args:
        file_ift            (str): name of intermediate flat file (".ift") including a reaction string
        verbose  (bool, optional): print more info regarding process status. Defaults to False.

    Returns:
        file_txt            (str): name of a modified intermediate flat file (".txt") including a reaction string
    """
    with open(file_ift) as f:
        for line in f.readlines():
            text = "".join(line)
        #text cleaning
            text2 = text.replace(",", ",\n").replace(":",":\t").replace("&#8652;","<->")
            file_txt = str(file_ift).replace(".ift", ".txt")
            g = open(file_txt, "w")
            g.write(text2)
            g.close()
            os.system("rm "+file_ift)

            if verbose:
                print("Info from file "+file_ift+" reordered")

    return file_txt

def get_string(file_txt, total_strings, verbose = False, remove_intermediate = True):
    """
    Helper function to extract reaction string and
    add to a reaction list in Reframed format.
    (Called in the "reframed_reaction_addition()" function).

    Args:
        file_txt                       (str): name of a modified intermediate flat file (".txt") including a reaction string
        total_strings                 (list): Python-list object including all reaction strings extracted from BiGG API processes
        verbose             (bool, optional): print more info regarding process status. Defaults to False.
        remove_intermediate (bool, optional): flag to remove intermediate ".ift" files. Defaults to True.
    """
    with open(file_txt) as asd:
        forbidden = ["+","<->","<--","-->"]
        for line in asd.readlines():
            line = line.strip()
            if line.startswith('"reaction_string"'):
                st = line.split("\t")[1]
                st1 = st.replace('"', "").replace("}", "").replace(',', '').strip()
                st2 = st1.split(" ")

                metabs = []
                for el in st2:
                    if not el in forbidden and not "." in el:
                        el_form = "M_"+el
                        metabs.append(el_form)
                    else:
                        metabs.append(el)

    #add Reframed name
    with open(file_txt) as asd:
        dsa = reversed(list(asd))
        for line in dsa:
            line = line.strip().replace("{", "")
            if line.startswith('"bigg_id"'):
                name = line.split("\t")[1]
                break
        right_name = "R_"+name.replace('"', "").replace(",", "").strip()+":"
        right_string = " ".join(metabs)
        reac_string = right_name+" "+right_string
        total_strings.append(reac_string)

        if verbose:
            print(str(file_txt)+" string PASSED") #debug

    #OPTIONAL - remove intermediate files
    if remove_intermediate:
        os.system("rm "+file_txt)

def get_string_error(file_txt, error_strings, verbose = False):
    """
    Helper function to extract reaction strings that lead to
    reaction addition errors (i.e. not found in BiGG).
    (Called in the "reframed_reaction_addition()" function).

    Args:
        file_txt                (str): name of a modified intermediate flat file (".txt") including a reaction string
        error_strings          (list): Python-list object including reaction strings that rose errors
        verbose      (bool, optional): print more info regarding process status. Defaults to False.
    """
    error_strings.append(file_txt)
    if verbose:
        print(str(file_txt)+" string NOT PASSED") #debug

def retry_genes(file_txt, verbose = False):
    """
    Helper function to try a different BiGG API request to
    extract reaction strings.
    (Called in the "reframed_reaction_addition()" function).

    Args:
        file_txt                (str): name of a modified intermediate flat file (".txt") including a reaction string
        verbose      (bool, optional): print more info regarding process status. Defaults to False.
    """
    cmd_line = "curl --silent 'http://bigg.ucsd.edu/api/v2/search?query=GENE&search_type=genes'"
    el = file_txt
    command2 = cmd_line.replace("GENE", el)
    file_ift2 = str(el)+".ift2"
    os.system(command2+" > "+file_ift2)
    if verbose:
        print(str(el)+" phase .ift2") #debug

    with open(file_ift2) as f:
        v = f.readline()
        v1 = v.split("}")[0].split("[")[1].replace("{","").split(", ")
        bigg_id = v1[0].replace('"', "").split(": ")[1]
        name = v1[1].replace('"', "").split(": ")[1]
        model_bigg_id = v1[3].replace('"', "").split(": ")[1]
        if not name == el:
            if verbose:
                print(str(el)+" string NOT PASSED AGAIN") #debug

            os.system("rm "+file_ift2)
            return
        else:
            cmd_line2 = "curl --silent 'http://bigg.ucsd.edu/api/v2/models/MODEL/genes/BIGG'"
            command3 = cmd_line2.replace("MODEL", model_bigg_id).replace("BIGG", bigg_id)
            file_ift3 = str(el)+".ift3"
            os.system(command3+" > "+file_ift3)
    os.system("rm "+file_ift2)

    with open(file_ift3) as f:
        v = f.readline().replace(",","\n")
        new_single_reac = v.split('reactions": [{"bigg_id": "')[1].split('"\n')[0]
    new_ift = curl_bigg_reaction(new_single_reac)
    new_txt = api_file_reorder2(new_ift, verbose = verbose)
    try:
        get_string(new_txt, total_strings, verbose = verbose)
        error_strings.remove(file_txt)
        os.system("rm "+file_ift3)
    except:
        return

def metab_change_names(model_new, verbose = False):
    """
    Helper function to find exact metabolites names
    and replacing the ones added via "get_string()".
    (Called in the "reframed_reaction_addition()" function).

    Args:
        model_new               (CBModel): Python-object copy of an existing ".xml" model, as generated via Reframed
        verbose          (bool, optional): print more info regarding process status. Defaults to False.
    """
    comm_list = []
    for metab in model_new.metabolites.values():
        metab_mod = str(metab).replace("M_", "", 1).replace("_c", "").replace("_e", "").replace("_p", "").replace("_m", "").replace("_n", "").replace("_x", "").replace("_h", "").replace("_g", "")
        if metab_mod in old_new_names.keys():
            command = "model_new.metabolites."+str(metab)+".name "+"="+' "'+old_new_names[metab_mod]+'"'
            comm_list.append(command)
            exec(command)

    if verbose:
        print(str(len(comm_list))+" commands used for METABOLITES name updates") #debug
    return

def check_name_changes(model_new):
    """
    Helper function to state the status of reactions included in the gap-filled model,
    i.e. wheter or not they've been renamed according to their proper name.
    (Called in the "reframed_reaction_addition()" function).

    Args:
        model_new   (CBModel): Python-object copy of an existing ".xml" model, as generated via Reframed
    """
    unchanged=[]
    for x in model_new.metabolites.values():
        if str(x).startswith("M_"):
            unchanged.append(x)

    if len(unchanged) == 0:
        result = "All METABOLITES names were updated"
    if len(unchanged) != 0:
        result = "ATTENTION Not all METABOLITES names were updated"
    print(result)

def reac_change_names(model_new, verbose = False):
    """
    Helper function to find the exact reaction names
    and replace the ones added in the new model via "get_string()".
    (Called in the "reframed_reaction_addition()" function).

    Args:
        model_new           (CBModel): Python-object copy of an existing ".xml" model, as generated via Reframed
        verbose     (bool, optional): print more info regarding process status. Defaults to False.
    """
    comm_list = []
    for reaz in model_new.reactions.values():
        reaz_mod = str(reaz.name).replace("R_", "", 1)
        if reaz_mod in old_new_names_R.keys():
            command = "model_new.reactions."+str(reaz.name)+".name "+"="+' "'+old_new_names_R[reaz_mod]+'"'
            try:
                exec(command)
                comm_list.append(command)
            except:
                command2 = "model_new.reactions."+str("R_"+reaz.name)+".name "+"="+' "'+old_new_names_R[reaz_mod]+'"'
                try:
                    exec(command2)
                    comm_list.append(command2)
                except:
                    pass
    if verbose:
        print(str(len(comm_list))+" commands used for REACTIONS name updates") #debug
    return

def old_new_names_dict(file):
    """
    Helper function to load info regarding correct names
    for BiGG metabolites.

    Args:
        file            (str): input ".tsv" file path

    Returns:
        old_new_names  (dict): Python-dictionary object:
                                keys: "uncorrected metabolite names"
                                values: "corrected metabolite names"
    """
    old_new_names = {}
    with open(file) as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            try:
                old = str(line[0])
                new = str(line[1])
                old_new_names.update({old:new})
            except:
                continue
    return old_new_names

def old_new_names_reac_dict(file):
    """
    Helper function to load info regarding correct names
    for BiGG reactions.

    Args:
        file                (str): input ".tsv" file path

    Returns:
        old_new_names_R    (dict): Python-dictionary object:
                                    keys: "uncorrected reaction names"
                                    values: "corrected reaction names"
    """
    old_new_names_R = {}
    with open(file) as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            try:
                old = str(line[0])
                new = str(line[1])
                old_new_names_R.update({old:new})
            except:
                continue
    return old_new_names_R

def fixed_modules_of_interest(dir_base, fixed_module_file):
    """
    Helper function to point out MODULES from a fixed list as in
    the command-line utilized instructions (".instruction" file).

    Args:
        dir_base            (str): folder path in which "kemet.py" was executed
        fixed_module_file   (str): file name of a ".instruction" file with indication of KEGG Modules of interest

    Returns:
        MODofinterest      (list): Python-list object including Modules of interest
    """
    MODofinterest = []

    os.chdir(dir_base)
    with open(fixed_module_file) as f:
        for line in f.readlines():
            MOD = line.strip()
            MODofinterest.append(MOD)
    return MODofinterest

def onbm_modules_of_interest(fasta_id, oneBM_modules_dir):
    """
    Helper function to points out MODULES with 1 block missing as in
    the command-line utilized instructions (".instruction" file).

    Args:
        fasta_id            (str): identificative FASTA name for a given MAG/Genome
        oneBM_modules_dir   (str): output folder of MAG/Genome specific ".instruction" file with KEGG Modules of interest

    Returns:
        MODofinterest      (list): Python-list object including Modules of interest
    """
    MODofinterest = []

    os.chdir(oneBM_modules_dir)
    if fasta_id+"_"+fixed_module_file in os.listdir():
        with open(fasta_id+"_"+fixed_module_file) as f:
            for line in f.readlines():
                MOD = line.strip()
                MODofinterest.append(MOD)
    return MODofinterest

def KOs_with_HMM_hits(hmm_hits_dir, fastakohits):
    """
    Helper function to point out KOs from modules identified with HMM HITS.

    Args:
        hmm_hits_dir    (str): HMM hits folder path
        fastakohits     (str): output file from "nhmmer_significant_hits_corr()"

    Returns:
        KOhits         (list): Python-list object including KOs identified from HMM procedures
    """
    KOhits = []

    os.chdir(hmm_hits_dir)
    with open(fastakohits) as f:
        for line in f.readlines()[1:]:
            KO = line.strip().split("\t")[0]
            KOhits.append(KO)
    return KOhits

def Modules_KOhits_connection(KOhits, MODofinterest):
    """
    Helper function to connect Modules with KOs missing, via dictionary.

    Args:
        KOhits                  (list): output of "KOs_with_HMM_hits()"
        MODofinterest           (list): Python-list object including Modules of interest

    Returns:
        MODofinterestXKOhits    (dict): Python-dictionary object:
                                        keys: "Modules of interest"
                                        values: "KOs missing in the Module"
    """
    MODofinterestXKOhits = {}

    for KO in KOhits:
        if KO == "":
            continue
        Modules = searchKeggShort(list_all_mod, Modules_directory, KO)
        for module in Modules:
            if not module in MODofinterest:
                continue
        # if KO is present in one of the modules of interest
            if not module in MODofinterestXKOhits.keys():
                vett_KOmiss = []
                vett_KOmiss.append(KO)
                MODofinterestXKOhits.update({module:vett_KOmiss})
        # if more than one KO is identified in the modules of interest
            else:
                vett_KOmiss = MODofinterestXKOhits[module]
                vett_KOmiss.append(KO)
                MODofinterestXKOhits.update({module:vett_KOmiss})
    return MODofinterestXKOhits

def modules_flat_files_connection(list_all_mod, MODofinterestXKOhits):
    """
    Connect Modules' names with the Modules flatfiles, via dictionary.

    Args:
        list_all_mod            (list): output of "list_all_modules()" function
        MODofinterestXKOhits    (dict): Python-dictionary output of "Modules_KOhits_connection()"

    Returns:
        modulesXflat            (dict): Python-dictionary object:
                                        keys: "KEGG Module id" - format: [Mxxxxx]
                                        values: "KEGG Module name"
    """
    modulesXflat = {}

    for module_txt in list_all_mod:
        for module in MODofinterestXKOhits.keys():
            if module_txt.startswith(module):
                modulesXflat.update({module:module_txt})
    return modulesXflat

def total_R_from_KOhits(MODofinterestXKOhits, modulesXflat, Modules_directory):
    """
    Point out RN from KOs, scanning in Modules flatfiles.

    Args:
        MODofinterestXKOhits    (dict): Python-dictionary output of "Modules_KOhits_connection()"
        modulesXflat            (dict): Python-dictionary output of "modules_flat_files_connection()"
        Modules_directory        (str): KEGG Modules flat-files folder path

    Returns:
        Rtotali_KOhits          (list): Python-list object including KEGG REACTIONS (R) defined by every KO identified from input dicts
    """
    Rtotali_KOhits = []
    KOhits_without_reaction = []

    os.chdir(Modules_directory)
    for module, module_txt in modulesXflat.items():
        missingKO = str(MODofinterestXKOhits[module]).replace("[","").replace("]","").replace("'","").split(", ")
        KOreac_module = knum4reac_mod(module_txt)
        for KO in missingKO:
            KO_s = KO.replace("(", "").replace(")", "")
            KO_split = re.split("[+-]", KO_s)
            for single_KO in KO_split:
                try:
                    reac_missing = KOreac_module[single_KO].split(",")
                    for single_reac in reac_missing:
                        if not single_reac in Rtotali_KOhits:
                            Rtotali_KOhits.append(single_reac)
                except:
                    KOhits_without_reaction.append(single_KO)
    return Rtotali_KOhits

def log_bigg_nr(bigg_nonredundant, fasta_id, gapfill_report_directory):
    """
    Helper function to save a log file of the output - FORMAT: 1 BiGG per line.

    Args:
        bigg_nonredundant          (list): Python-list output of "bigg_nonredundant()" function
        fasta_id                    (str): identificative FASTA name for a given MAG/Genome
        gapfill_report_directory    (str): output folder path
    """
    os.chdir(gapfill_report_directory)
    f = open("bigg_log_"+fasta_id+".txt", "w")
    for single_bigg in bigg_nonredundant:
        f.write(single_bigg+"\n")
    f.close()
    print("COMPLETE BiGG logging "+fasta_id) #debug

def log_modelseed_nr(modelseed_nonredundant, fasta_id, gapfill_report_directory):
    """
    Helper function to save a log file of the output - FORMAT: 1 SEED per line.

    Args:
        modelseed_nonredundant     (list): Python-list output of "modelseed_nonredundant()" function
        fasta_id                    (str): identificative FASTA name for a given MAG/Genome
        gapfill_report_directory    (str): output folder path
    """
    os.chdir(gapfill_report_directory)
    f = open("seed_log_"+fasta_id+".txt", "w")
    for single_rxn in modelseed_nonredundant:
        f.write(single_rxn+"\n")
    f.close()
    print("COMPLETE ModelSEED logging "+fasta_id) #debug

def reframed_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = False):
    """
    Main function to add reactions linked to HMM-derived genomic evidence
    obtained from prior KEMET processes.
    Makes large use of BiGG API to recover info regarding reaction strings.
    Include said reaction strings in a copy of an existing genome-scale model
    through the use of Reframed package functions.

    Args:
        fasta_id                        (str): identificative FASTA name for a given MAG/Genome
        model_directory                 (str): existing CarveMe model folder path
        gapfill_report_directory        (str): report output folder path
        bigg_api                        (str): folder path in which to store files with reaction strings (via BiGG API)
        verbose              (bool, optional): print more info regarding process status. Defaults to False.
    """
    datetoday = str(datetime.now())[:10]
    # MODEL IO VIA REFRAMED
    os.chdir(model_directory)
    for file in sorted(os.listdir()):
        if file.endswith(".xml") and file.startswith(fasta_id):
            os.chdir(model_directory)
            model = load_cbmodel(file)

    # Task: saves reaction to add (from "bigg_nonredundant()") in a list
            os.chdir(gapfill_report_directory)
            bigg_gapfill = bigg_gapfill_absent_in_model("bigg_log_"+fasta_id+".txt", model)

    # Task: create a directory to store downloaded reaction from BiGG API
            os.chdir(bigg_api)
            new_modelgapfill_directory = model.id+"_gapfill_directory"
            if not new_modelgapfill_directory in os.listdir():
                os.mkdir(new_modelgapfill_directory)

            modelgapfill_directory = bigg_api+new_modelgapfill_directory

    # ACTUAL REACTION RESOURCE DOWNLOAD
            os.chdir(modelgapfill_directory)

            total_strings=[]
            error_strings=[]
            if __name__ == '__main__':
                with Pool(processes=6) as p: # POSSIBILITY: change number of processes/threads (default = 6)
                    p.map(curl_bigg_reaction, bigg_gapfill)
                    p.close()

            lista = os.listdir()
            for file in lista:
                if file.endswith(".ift"):
                    api_file_reorder(file, verbose = verbose)
            lista = sorted(os.listdir())
            for file in lista:
                if file.endswith(".txt"):
                    try:
                        get_string(file, total_strings, verbose = verbose)
                    except:
                        get_string_error(file, error_strings, verbose = verbose)
                        try:
                            retry_genes(file.replace(".txt",""), verbose = verbose)
                        except:
                            pass
    # Task: adds the reactions to a copy of the model, without those having any non-prokaryotic compartments
            model_new = model.copy()
            exo_comp_reacs = []
            added_reacs = []

            for reac in total_strings:
                reacs = reac.strip().split(" ")
                for element in reacs:
                    if element.endswith("_n") or element.endswith("_g") or element.endswith("_h") or element.endswith("_x") or element.endswith("_m") or element.endswith("_r") or element.endswith("_v"):
                        exo_comp_reacs.append(reac)
                        break

            for reac in total_strings:
                if reac in exo_comp_reacs:
                    continue
                model_new.add_reaction_from_str(reac)
                added_reacs.append(reac)

                if verbose:
                    print("ADDED "+reac) #debug

    # Task:  saves a LOGfile with the reactions ACTUALLY ADDED to the model
            os.chdir(gapfill_report_directory)
            f = open(""+model.id+"_added_reactions.txt", "w")
            for reac in added_reacs:
                f.write(reac+"\n")
            f.close()

    # Task:  checks namespace quality for metabolites and saves the new gap-filled model
            if verbose:
                metab_change_names(model_new, verbose = verbose)
                reac_change_names(model_new, verbose = verbose)
                check_name_changes(model_new)

            os.chdir(gapfilled_model_directory)
            save_cbmodel(model_new, model_new.id+"_KEGGadd_"+datetoday+".xml", flavor = 'bigg')

def recap_addition(fasta_id, gapfill_report_directory, old_new_names_R):
    """
    Stores RECAP informations for reaction addition
    (via "existing" GSMM procedures, specified in command-line args) as performed in each model.
    This function generates a final tabular output of the GSMM process.

    Args:
        fasta_id                    (str): identificative FASTA name for a given MAG/Genome
        gapfill_report_directory    (str): output folder path
        old_new_names_R            (dict): Python-dictionary output of "old_new_names_reac_dict()" function
    """
    datetoday = str(datetime.datetime.now())[:10]

    os.chdir(gapfill_report_directory)
    f = open("recap_gapfill_"+datetoday+".tsv", "a")
    header = "MAG\tSTRING\tNAME\n"
    if f.tell() == 0:
        f.write(header)

    for file in sorted(os.listdir()):
        if file.endswith("_added_reactions.txt") and fasta_id in file:
            with open(file) as g:
                MAG = file[:-20]
                for line in g.readlines():
                    STRING = line.strip()
                    reac_id = STRING.split(": ")[0].replace("R_", "")
                    NAME = old_new_names_R[reac_id]
                    if not STRING == "":
                        f.write(MAG+"\t"+STRING+"\t"+NAME+"\n")
    f.close()

########################################################################################
if __name__ == "__main__":
    import os
    import argparse
    from datetime import datetime

    run_start = str(datetime.now().strftime("%Y-%m-%d"))
    LOGflag = False

    ###############
    # directories #
    ###############
    dir_base = os.getcwd() #script folder

    Modules_directory = dir_base+"/KEGG_MODULES/"
    kkfiles_directory = Modules_directory+"/kk_files/"
    KAnnotation_directory = dir_base+"/KEGG_annotations/"
    output_directory = dir_base

    taxa_dir = dir_base+"/taxonomies/"
    dir_base_KO = dir_base+"/Knumber_ntsequences/"
    msa_dir = dir_base+"/multiple_fasta/"
    hmm_dir = dir_base+"/HMM/"
    hmm_hits_dir = dir_base+"/HMM_HITS/"
    dir_genomes = dir_base+"/genomes/"
    oneBM_modules_dir = dir_base+"/oneBM_modules/"
    gapfill_report_directory = dir_base+"/report_gapfill/"
    bigg_api = dir_base+"/biggapi_download/"
    DB_directory = dir_base+"/DB/"
    model_directory = dir_base+"/models/"
    gapfilled_model_directory = dir_base+"/models_gapfilled/"
    de_novo_model_directory = dir_base+"/de_novo_models/"

    #####################
    # INSTRUCTION_FILES #
    #####################
    instruction_file = "genomes.instruction"
    fixed_module_file = "module_file.instruction"
    fixed_ko_file = "ko_file.instruction"

    #################
    # BASE COMMANDS #
    #################
    base_com_KEGGget = _base_com_KEGGget
    base_com_mafft = _base_com_mafft
    base_com_hmmbuild = _base_com_hmmbuild
    base_com_nhmmer = _base_com_nhmmer

    #############################################

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description =
    '''
    KEMET pipeline:
    1) Evaluate KEGG Modules Completeness for given genomes.
    2) HMM-based check for ortholog genes (KO) of interest after KEGG Module Completeness evaluation.
    3) Genome-scale model gapfill with nucleotidic HMM-derived evidence, for KOs of interest.
    ''')

    parser.add_argument('FASTA_file',
                        help='''Genome/MAG FASTA file as indicated in the "genomes.instruction" -
    points to files (in "KEGG_annotations") comprising KO annotations, associated with each gene.''')
    parser.add_argument('-a','--annotation_format', required=True,
                        choices=_ktest_formats,
                        help='''Format of KO_list.
    eggnog: 1 gene | many possible annotations;
    kaas: 1 gene | 1 annotation at most;
    kofamkoala: 1 gene | many possible annotations''')
    parser.add_argument('--update_taxonomy_codes', action ="store_true",
                        help='''Update taxonomy filter codes - WHEN TO USE: after downloading a new BRITE taxonomy with "setup.py".''')
    parser.add_argument('-I','--path_input',
                        help='''Absolute path to input file(s) FOLDER.''', default = KAnnotation_directory)
    parser.add_argument('-k','--as_kegg', action ="store_true",
                        help='''Return KEGG-Mapper output for the Module Completeness evaluation.''')
    parser.add_argument('--skip_hmm', action ="store_true",
                        help='''Skip HMM-driven search for KOs & stop after KEGG Modules Completeness evaluation.''')
    parser.add_argument('--hmm_mode',
                        choices=_hmm_modes,
                        help='''Choose the subset of KOs of interest for HMM-based check.
    By default, the KOs already present in the functional annotation are not checked further.

    onebm: search for KOs from KEGG Modules missing 1 block;
    modules: search for KOs from the KEGG Modules indicated in the "module_file.instruction" file, 1 per line
        (e.g. Mxxxxx);
    kos: search for KOs indicated in the "ko_file.instruction" file, 1 per line
        (e.g. Kxxxxx)
                        ''')
    parser.add_argument('--threshold_value', default = _def_thr,
                        help='''Define a threshold for the corrected score resulting from HMM-hits, which is indicative of good quality.''')
    parser.add_argument('--skip_nt_download', action ="store_true",
                        help='''Skip downloading KEGG KOs nt sequences.''')
    parser.add_argument('--skip_msa_and_hmmbuild', action ="store_true",
                        help='''Skip MAFFT and HMMER hmmbuild commands.''')
    parser.add_argument('--retry_nhmmer', action ="store_true",
                        help='''Move HMM-files and re-run nHMMER command.''')

    parser.add_argument('--skip_gsmm', action ="store_true",
                        help='''Skip GSMM operations, gapfill or de-novo model creation, & stop after HMM-driven search for KOs.''')
    parser.add_argument('--gsmm_mode',
                        choices=_gapfill_modes,
                        help='''Choose the methods of GSMM operation.
    (This method won't be performed if "--hmm_mode kos" was chosen)
    existing: use pre-existing CarveMe GSMM to add reactions content connected to HMM-derived KOs;
    denovo: generate a new CarveMe GSMM, performing gene prediction and adding HMM-derived hits from the chosen HMM-mode.
                        ''')
    parser.add_argument('-O','--path_output',
                        help='''Absolute path to ouput file(s) FOLDER.''', default = dir_base)
    parser.add_argument('-v','--verbose', action ="store_true",
                        help='''Print more informations - for debug and progress.''')
    parser.add_argument('-q','--quiet', action ="store_true",
                        help='''Silence soft-errors (for MAFFT and HMMER commands).''')
    parser.add_argument('--log', action ="store_true",
                    help='''Store KEMET commands and progress during the execution in a log file.''')
    args = parser.parse_args()

#### SET NEW IO FOLDERS
    KAnnotation_directory = args.path_input
    output_directory = args.path_output
    report_txt_directory = output_directory+"/reports_txt/"
    report_tsv_directory = output_directory+"/reports_tsv/"
    ktests_directory = output_directory+"/ktests/"
    klists_directory = output_directory+"/klists/"
    if args.log:
        LOGflag = True
        import logging
        logging.basicConfig(filename=run_start+'_KEMET_execution.log',
        level=logging.INFO, format='%(asctime)s %(message)s')
        KEMET_args = str(args).replace('Namespace(','').replace(')','')
        logging.info('KEMET arguments: '+KEMET_args)

#### KMC - PRODUCE KTEST FILE
    os.chdir(KAnnotation_directory)
    file_name = str(args.FASTA_file).rsplit("/",1)[-1].replace(".fasta","").replace(".fna","").replace(".fa","") # from path indication of a contig file maintain file_name
    if LOGflag:
        logging.info('+++\tSTART '+file_name)
        logging.info('+++START KEGG Modules completeness')
    for file in os.listdir():
        if file.startswith(file_name):
            if args.annotation_format == "kaas":
                if args.verbose:
                    print("converting kaas-like file "+file)
                ktest, KOs = KAASXktest(file, file.rsplit(".",1)[0]+".ktest", KAnnotation_directory, ktests_directory)
            elif args.annotation_format == "eggnog":
                if args.verbose:
                    print("converting eggnog file "+file)
                ktest, KOs = eggnogXktest(file, file.rsplit(".",2)[0]+".ktest", KAnnotation_directory, ktests_directory)
            elif args.annotation_format == "kofamkoala":
                if args.verbose:
                    print("converting kofamkoala file "+file)
                ktest, KOs = kofamXktest(file, file.rsplit(".",1)[0]+".ktest", KAnnotation_directory, ktests_directory)
            else:
                raise('''
                The only accepted formats for the --annotation_format method are: {}".format(_ktest_formats)
                ''')

    os.chdir(ktests_directory)
    if ktest in sorted(os.listdir()):
        ko_list = create_KO_list(ktest, ktests_directory)
        os.chdir(kkfiles_directory)
        for file in sorted(os.listdir()):
            if file.endswith(".kk"):
                testcompleteness(ko_list, file, kkfiles_directory, report_txt_directory, "reportKMC_"+ktest[:-6]+".txt")
                if args.as_kegg:
                    testcompleteness_tsv(ko_list, file, kkfiles_directory, report_tsv_directory, "reportKMC_"+ktest[:-6]+".tsv", as_kegg = True)
                else:
                    testcompleteness_tsv(ko_list, file, kkfiles_directory, report_tsv_directory, "reportKMC_"+ktest[:-6]+".tsv")
        if LOGflag:
            logging.info('COMPLETE KEGG Modules completeness')

    if args.skip_hmm:
        raise SystemExit('''
        You have chosen not to perform HMM-driven ortholog search & later analyses.
        If you didn't want to stop here, remove the --skip_hmm argument from the command line.
        ''')

#### HMM - VERBOSITY SETTINGS
    if args.verbose:
        base_com_mafft = base_com_mafft.replace("--quiet","")
        base_com_hmmbuild = base_com_hmmbuild.replace(" > /dev/null","")
        base_com_nhmmer = base_com_nhmmer.replace(" > /dev/null","")
    if args.quiet:
        base_com_mafft = base_com_mafft+" 2>/dev/null"
        base_com_hmmbuild = base_com_hmmbuild+" 2>&1"
        base_com_nhmmer = base_com_nhmmer+" 2>&1"

#### HMM - READ AND WRITE INSTRUCTIONS
    os.chdir(dir_base)
    FASTA = str(args.FASTA_file).rsplit("/",1)[-1]
    with open(instruction_file) as f:
        for line in f.readlines()[1:]:
            if line.startswith(FASTA):
                line = line.strip().split("\t")
                if len(line) == 3:
                    metabolic_universe = line[2]
                for file in os.listdir(dir_genomes):
                    if file.startswith(FASTA):
                        fasta_genome = dir_genomes+FASTA

                fasta_id = FASTA.replace(".fasta","").replace(".fna","").replace(".fa","")
                klist_file = fasta_id+".klist"
                taxonomy = line[1]
                taxa_file = taxonomy+".keg"
                if args.hmm_mode == "modules":
                    os.chdir(dir_base)
                    tuple_modules = create_tuple_modules(fixed_module_file)
                    write_KOs_from_modules(fasta_id, tuple_modules, report_txt_directory, klists_directory)
                if args.hmm_mode == "onebm":
                    tuple_modules = create_tuple_modules_1BM(fasta_id, fixed_module_file, oneBM_modules_dir, report_tsv_directory)
                    write_KOs_from_modules(fasta_id, tuple_modules, report_txt_directory, klists_directory)
                if args.hmm_mode == "kos":
                    write_KOs_from_fixed_list(fasta_id, fixed_ko_file, ktests_directory, klists_directory)

                dir_KO = dir_base_KO+taxonomy.replace(" ","_")+"/"                           # KO folder for taxonomy,    save time!
                msa_dir_comm = msa_dir+fasta_id+"/"                                          # MSA folder for MAG/Genome, more ordered!
                hmm_dir_comm = hmm_dir+fasta_id+"/"                                          # HMM folder for MAG/Genome, more ordered!
                CORR_THRESHOLD = float(args.threshold_value)

#### HMM - OPERATE SINGLE FUNCTIONS
                print(_timeinfo()+"+++START HMM operations "+fasta_id)
                if LOGflag:
                    logging.info('+++START HMM operations '+fasta_id)
                if args.update_taxonomy_codes:
                    taxa_allow = taxonomy_filter(taxonomy, dir_base, taxa_file, taxa_dir, update = True)
                else:
                    os.chdir(taxa_dir)
                    if not taxa_file in os.listdir():
                        taxa_allow = taxonomy_filter(taxonomy, dir_base, taxa_file, taxa_dir, update = True)
                    else:
                        taxa_allow = taxonomy_filter(taxonomy, dir_base, taxa_file, taxa_dir)

                if not args.skip_nt_download:
                    if LOGflag:
                        logging.info('START download nucleotidic sequences')
                    download_ntseq_of_KO(klist_file, dir_base_KO, dir_KO, klists_directory, taxa_dir, taxa_file, base_com_KEGGget)
                    if LOGflag:
                        logging.info('COMPLETE download nucleotidic sequences')
                if args.retry_nhmmer:
                    # POSSIBILITY: after a whole KEMET run, to try other parameters
                    movebackHMM(hmm_dir_comm, msa_dir_comm)
                if not args.skip_msa_and_hmmbuild:
                    if LOGflag:
                        logging.info('START sequences filtering and allignment')
                    filter_and_allign(taxa_dir, taxa_file, fasta_id, klist_file, klists_directory, msa_dir, dir_KO)
                    if LOGflag:
                        logging.info('COMPLETE Filter and allign')
                    if LOGflag:
                        MSA_and_HMM(msa_dir_comm, base_com_mafft, base_com_hmmbuild, log=True)
                    else:
                        MSA_and_HMM(msa_dir_comm, base_com_mafft, base_com_hmmbuild, log=False)
                nhmmer_for_genome(fasta_genome, msa_dir_comm, base_com_nhmmer)
                if LOGflag:
                    logging.info('COMPLETE nhmmer')
                move_HMM_and_clean(hmm_dir_comm, msa_dir_comm)

#### HMM - FIRST REPORT FILE
                nhmmer_significant_hits_corr(fasta_id, hmm_dir_comm, corr_threshold = CORR_THRESHOLD)
                if LOGflag:
                    logging.info('COMPLETE nhmmer significant hits')
                HMM_hits_dict = HMM_hits_sequences(hmm_hits_dir, dir_genomes)
                HMM_hits_TRANSLATED_dict = HMM_hits_translated_sequences(HMM_hits_dict)
                HMM_hits_longestTRANSLATED_dict = HMM_hits_longest_translated_sequences(HMM_hits_dict, HMM_hits_TRANSLATED_dict)

#### HMM - TOTAL REPORT FILE
                recap_hits_corr(fasta_id, hmm_hits_dir, HMM_hits_dict, HMM_hits_longestTRANSLATED_dict, run_start)

                if args.skip_gsmm:
                    print('''
                    You have chosen not to perform GSMM gapfilling.
                    If you didn't want to stop here, remove the --skip_gsmm argument from the kemet.py command line.
                    ''')
                else:
                    if args.hmm_mode == "kos" and args.gsmm_mode == "existing":
                        print('''
                        The GSMM operations chosen are not compatible with the HMM mode, therefore no gapfilling would be performed.
                        ''')
                    else:
#### GSMM - OPERATE SINGLE FUNCTIONS
                        if args.gsmm_mode == "denovo" or args.gsmm_mode == "existing":
                            if LOGflag:
                                logging.info('+++START GSMM operations '+fasta_id)

                            fastakohits = FASTA+"_HMM_hits.txt"
                            current_run = "KEMET_run_"+run_start

                            if args.gsmm_mode == "denovo":
                                if LOGflag:
                                    build_de_novo_GSMM(FASTA, fasta_genome, de_novo_model_directory, current_run, log=True)
                                else:
                                    build_de_novo_GSMM(FASTA, fasta_genome, de_novo_model_directory, current_run, log=False)

                            elif args.gsmm_mode == "existing":
                                list_all_mod = list_all_modules(Modules_directory)
                                os.chdir(DB_directory)
                                DB_KEGG_RN = KEGG_BiGG_SEED_RN_dict("reactions_DB.tsv", DB_directory, ontology="BiGG")
                                old_new_names = old_new_names_dict("metabolites_names_from_id_bigg.tsv")
                                old_new_names_R = old_new_names_reac_dict("reactions_names_from_id_bigg.tsv")

                                if args.hmm_mode == "modules":
                                    MODofinterest = fixed_modules_of_interest(dir_base, fixed_module_file)
                                if args.hmm_mode == "onebm":
                                    MODofinterest = onbm_modules_of_interest(fasta_id, oneBM_modules_dir)
                                KOhits = KOs_with_HMM_hits(hmm_hits_dir, fastakohits)
                                MODofinterestXKOhits = Modules_KOhits_connection(KOhits, MODofinterest)
                                modulesXflat = modules_flat_files_connection(list_all_mod, MODofinterestXKOhits)

                                Rtotali_KOhits = total_R_from_KOhits(MODofinterestXKOhits, modulesXflat, Modules_directory)

                                KEGG_R_to_add = keggR_in_DB(Rtotali_KOhits, DB_KEGG_RN)
#### GSMM - REACTION ADDITION + RECAP
                                bigg_nonredundant = bigg_nonredundant(KEGG_R_to_add, DB_KEGG_RN)
                                log_bigg_nr(bigg_nonredundant, fasta_id, gapfill_report_directory)
                                if args.verbose:
                                    reframed_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = True)
                                else:
                                    reframed_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = False)
                                recap_addition(fasta_id, gapfill_report_directory, old_new_names_R)

                print(_timeinfo()+"END "+fasta_id+"\n")
                if LOGflag:
                    logging.info('END '+fasta_id+'\n')