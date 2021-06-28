#!/usr/bin/env python
# coding: utf-8

import os
import re

###############
# extra specs #
###############
_ktest_formats = ["eggnog", "kaas", "kofamkoala"] #TODO: add other formats, BlastKOALA-GhostKOALA
_output_formats = ["txt", "tsv", "txt+tsv"]
###############

def eggnogXktest(eggnog_file, converted_output, KAnnotation_directory, ktests_directory):
    """ Starting from eggNOG output (1 gene - many annotations), keep only KO """
    os.chdir(KAnnotation_directory)
    KOs = {}

    with open(eggnog_file) as g:
        headers = g.readlines()[3].strip().split("\t")
        koslice = 0
        for field in headers:
            if not field == "KEGG_ko":
                koslice += 1
            else:
                break
        g.seek(0)
        for line in g.readlines()[4:-3]: # skip info lines w/o genes
            fasta_id = line.strip().split("\t")[0]
            egg_kos = line.strip().split("\t")[koslice].replace("ko:","")
            #print(fasta_id) # debug
            #print(egg_kos) # debug
            if egg_kos != "":
                egg_kos_hits = egg_kos.split(",")
                for ko in egg_kos_hits:

                    if not ko in KOs:
                        #print(ko) #debug
                        KOs[ko] = 1
                    else:
                        KOs[ko] += 1

        # POSSIBILITY: for each gene, correcting per diff. ortholog hits - if more KOs -> fraction of KO cound

                    #if not ko in KOs:
                        #KOs[ko] = round(1/len(egg_kos_hits), 2) # correction
                    #else:
                        #KOs[ko] += round(1/len(egg_kos_hits), 2) # correction
            else:
                pass
    #print(str(len(KOs))) # debug

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
    """ Starting from KAAS output (1 gene - 1 KO), keep only KO """
    os.chdir(KAnnotation_directory)
    KOs = []
    with open(file_kaas) as f:
        for line in f.readlines():
            line_s = line.strip().split("\t")
            if len(line_s) == 2:
                if not line_s[1] in KOs: ## non-redundant KO addition
                    # POSSIBILITY: KO for each gene, storing them all
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
    """ Starting from KofamKOALA output (1 gene - many annotations), keep only KO """
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
            #print(fasta_id) # debug
            #print(kofam_ko) # debug
            if kofam_ko != "":
                if not kofam_ko in KOs:
                    #print(kofam_ko) #debug
                    KOs[kofam_ko] = 1
                else:
                    KOs[kofam_ko] += 1

        # POSSIBILITY: for each gene, correcting per diff. ortholog hits - if more KOs -> fraction of KO cound
            # ADD CODE SIMILAR TO eggnogXktest
            else:
                pass

    #print(str(len(KOs))) # debug
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
    """ From KOs file (.ktest file), return a Python list """
    os.chdir(ktests_directory)
    ko_list = []
    with open(file_ko_list) as f:
        for line in f.readlines():
            line_s = line.strip()
            ko_list.append(line_s)
    return ko_list

def testcompleteness(ko_list, kk_file, kkfiles_directory, report_txt_directory, file_output = "report.txt", cutoff = 0):
    '''
    ko_list: non-redundant ko list, produced from KAAS-like/eggNOG/other annotators - 1 KO x line 
    kk_file: module file ".kk", with KO indication for every block - includes COMPLEX and OPTIONAL KOs
    kkfiles_directory: directory of ".kk" files - to be scanned
    report_txt_directory: output directory for .txt flatfile KMC report
    file_output: output file name
    cutoff: optional parameter, in order to set which is the least % of the modules to be included in the output
    '''
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
            #print(complexes) # debug
# search optionals
        if v_optional != -1:
            o_list = v[v_optional+1].replace("\n", "").replace("\t", "").split(", ")
            for el in o_list:
                optional.append(el)
            ko_list_optional = [el for el in ko_list]
            for el in optional:
                ko_list_optional.append(el)
            #print(optional) # debug

        for line in v[1:end]:
            count_lines += 1
            #print(line) #debug
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
                        #print(k_singlecomplex) # debug
                        if all(el in ko_line for el in k_singlecomplex): # if EACH complex-part in line
                            if all(el in ko_list_optional for el in k_singlecomplex):  # if element from KOlist+optional
                                check = 1
                            else:
                                continue
                    else:
                        for element in ko_line:
                            if not element in str(complexes):
                                if element in ko_list:
                                    #print(element) # debug
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
                    #print(control) # debug
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
                        #print(control) # debug
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
                missing_blocks=str(present)+"__"+str(total)
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

    if percentage_round >= cutoff: # optional parameter
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
    #return report # debug

def testcompleteness_tsv(ko_list, kk_file, kkfiles_directory, report_tsv_directory, file_report_tsv = "report.tsv", as_kegg = False, cutoff = 0):
    '''
    ko_list: non-redundant ko list, produced from KAAS-like/eggNOG/other annotators - 1 KO x line 
    kk_file: module file ".kk", with KO indication for every block - includes COMPLEX and OPTIONAL KOs
    kkfiles_directory: directory of ".kk" files - to be scanned
    report_tsv_directory: output directory for .tsv tabular KMC report
    file_output: output file name
    cutoff: optional parameter, in order to set which is the least % of the modules to be included in the output
    '''
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
            #print(complexes) # debug

    # list optionals
        if v_optional != -1:
            o_list = v[v_optional+1].replace("\n", "").replace("\t", "").split(", ")
            for el in o_list:
                optional.append(el)
            ko_list_optional = [el for el in ko_list]
            for el in optional:
                ko_list_optional.append(el)
            #print(optional) # debug

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
                    #print(KO) # debug
                    Kpresent.append(KO)
            else:
                Kmissing.append(KO)

    # CHECKS for each line in .kk file: KOs, complexes, optionals and list-indication
        for line in v[1:end]:
            count_lines += 1
            #print(line.strip()) #debug
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
                        if all(el in ko_line for el in k_singlecomplex): # if EACH complex-part in line
                            if all(el in ko_list_optional for el in k_singlecomplex):  # if element from KOlist+optional

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

        #percentage_2digits_tsv = percentage_2digits
        #if percentage_2digits == 100:
        #    completeness = "COMPLETE"
        #else:
        #    completeness = "INCOMPLETE"

        #print(str(count_lines)+" "+line.strip()+" "+missing_blocks)
    # REPORT INFOS

        ## POSSIBILITY: write 1-2 BLOCK(S) MISSING only for >3 blocks Modules, as KEGG output
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

        #report.insert(1, "%\t"+str(percentage_2digits)+"\t"+str(present)+"__"+str(present+missing)+"\t"+completeness+"\n")
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
            #print(el) # debug
            if type(el) == str:
                h.write(el+"\t")
            if type(el) == list:
                knums = ",".join(el)
                h.write(knums+"\t")
        h.write("\n")
        h.close()

    #return report # debug

if __name__ == "__main__":
    import os
    import argparse

    ###############
    # directories #
    ###############
    dir_base = os.getcwd()

    Modules_directory = dir_base+"/KEGG_MODULES/"
    kkfiles_directory = Modules_directory+"/kk_files/"
    KAnnotation_directory = dir_base+"/KEGG_mappings/"
    output_directory = dir_base

    ###############

    parser = argparse.ArgumentParser(description="Evaluate KEGG Modules Completeness for given genomes.")

    parser.add_argument('KOfile',
                        help='''File comprising KO annotations, associated with each gene.''')
    parser.add_argument('-a','--annotation_format', required=True,
                        choices=_ktest_formats,
                        help='''Format of KO_list.
                        eggnog: 1 gene | many possible annotations;
                        kaas: 1 gene | 1 annotation at most.''')
    parser.add_argument('-I','--path_input',
                        help='''Absolute path to input file(s) FOLDER.''', default = KAnnotation_directory)
    parser.add_argument('-o','--output', default = "tsv",
                        choices=_output_formats,
                        help='''Output format for KMC summary.
                        txt: more level-detailed, worse in recap;
                        tsv: best at recap, easily parsable for downstream analysis;
                        txt+tsv: both of the above.''')
    parser.add_argument('-k','--as_kegg', action ="store_true",
                        help='''Return KEGG Mapper output for the Module completeness.''')
    parser.add_argument('-O','--path_output',
                        help='''Absolute path to ouput file(s) FOLDER.''', default = dir_base)
    parser.add_argument('-v','--verbose', action ="store_true",
                        help='''Print more informations - for debug and progress.''')
    args = parser.parse_args()
    #print(args) #debug

#### SET NEW IO FOLDERS

    KAnnotation_directory = args.path_input
    output_directory = args.path_output
    report_txt_directory = output_directory+"/reports_txt/"
    report_tsv_directory = output_directory+"/reports_tsv/"
    ktests_directory = output_directory+"/ktests/"

#### PRODUCE KTEST FILE

    os.chdir(KAnnotation_directory)
    file = args.KOfile
    if file in os.listdir():
        if args.annotation_format == "kaas":
            if args.verbose:
                print("converting kaas-like file "+file)
            #POSSIBILITY: adding "ko_list_" in intermediate file - needs slices afterwards
            ktest, KOs = KAASXktest(file, file.rsplit(".",1)+".ktest", KAnnotation_directory, ktests_directory)
        elif args.annotation_format == "eggnog":
            if args.verbose:
                print("converting eggnog file "+file)
            ktest, KOs = eggnogXktest(file, file.rsplit(".",2)+".ktest", KAnnotation_directory, ktests_directory)
        elif args.annotation_format == "kofamkoala":
            if args.verbose:
                print("converting kofamkoala file "+file)
            ktest, KOs = kofamXktest(file, file.rsplit(".",1)+".ktest", KAnnotation_directory, ktests_directory)
        else:
            raise("The only accepted formats for the --annotation_format method are the one made by: {}".format(_ktest_formats))

    os.chdir(ktests_directory)
    if ktest in sorted(os.listdir()):
        ko_list = create_KO_list(ktest, ktests_directory)
        os.chdir(kkfiles_directory)
        for file in sorted(os.listdir()):
            if file.endswith(".kk"):
#### TESTCOMPLETENESS .txt
                if args.output == "txt":
                    testcompleteness(ko_list, file, kkfiles_directory, report_txt_directory, "reportKMC_"+ktest[:-6]+".txt")
#### TESTCOMPLETENESS .tsv
                if args.output == "tsv":
                    if args.as_kegg:
                        testcompleteness_tsv(ko_list, file, kkfiles_directory, report_txt_directory, "reportKMC_"+ktest[:-6]+".tsv", as_kegg = True)
                    testcompleteness_tsv(ko_list, file, kkfiles_directory, report_tsv_directory, "reportKMC_"+ktest[:-6]+".tsv")
#### TESTCOMPLETENESS BOTH .txt AND .tsv
                if args.output == "txt+tsv":
                    testcompleteness(ko_list, file, kkfiles_directory, report_txt_directory, "reportKMC_"+ktest[:-6]+".txt")
                    if args.as_kegg:
                        testcompleteness_tsv(ko_list, file, kkfiles_directory, report_txt_directory, "reportKMC_"+ktest[:-6]+".tsv", as_kegg = True)
                    testcompleteness_tsv(ko_list, file, kkfiles_directory, report_tsv_directory, "reportKMC_"+ktest[:-6]+".tsv")

#### NUMBER OF COMPLETE MODULES PER MODEL
    if args.verbose and args.output in ("txt","txt+tsv"):
        os.chdir(report_txt_directory)
        filereport = "reportKMC_"+ktest[:-6]+".txt"
        if filereport in sorted(os.listdir()):
            with open(filereport) as f:
                n = 0
                for line in f.readlines():
                    if "\tCOMPLETE" in line:
                        n += 1
                else:
                    print(filereport+"\t"+str(n))
