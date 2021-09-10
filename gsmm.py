#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import os
import re
import argparse
import datetime
from multiprocessing import Process
from multiprocessing import Pool
#from cobra.io import read_sbml_model #TODO: implement
import reframed
from reframed import load_cbmodel
from reframed import save_cbmodel

def build_de_novo_GSMM(FASTA, fasta_genome, de_novo_model_directory, hmm_hits_dir):
    '''
    Generate de-novo GSMM from Prodigal gene calling, adding HMM-hits derived from hmm.py
    ------------------------------------------
    INPUTS:  FASTA              - name of MAG of interest
             fasta_genome       - name and path of MAG of interest
             <file_recap_> file - HMM-hits recap file
    OUTPUTS: gene prediction, DIAMOND intermediates, GSMM of MAG of interest (.xml)
    '''
    HMM_HITS = {}
    os.chdir(hmm_hits_dir)
    for file in os.listdir():
        if not file.startswith("file_recap_") and file.endswith(".tsv"): # open the recap from hmmm.py
            continue
        with open(file) as f:
            for line in f.readlines()[1:]:
                line = line.strip().split("\t")
                MAG = line[0]
                KO = line[1]
                xseq = line[12]
                if not MAG == FASTA:
                    continue
                HMM_HITS.update({">"+KO:xseq})

    # launch prodigal gene calling
    os.chdir(de_novo_model_directory)
    if not "proteins" in os.listdir():
        os.mkdir("proteins")
    prodigal_cmd = "prodigal -i fasta_genome -a ./proteins/FASTA.faa -q -m"
    prodigal_cmd_mod = prodigal_cmd.replace("fasta_genome",fasta_genome).replace("FASTA",FASTA)
    os.system(prodigal_cmd_mod)

    # ADD HMM-hits
    os.chdir("./proteins")
    f = open(FASTA+".faa", "a")
    for ko, xseq in HMM_HITS.items():
        f.write(ko+"\n")
        f.write(xseq+"\n")
    f.close()
    os.chdir("..")

    # GENERATE MODEL
    if not "dmnd_intermediates" in os.listdir():
        os.mkdir("dmnd_intermediates")
    carveme_cmd = "carve ./proteins/FASTA.faa --fbc2 -u universe -o FASTA.xml"
    carveme_cmd_mod = carveme_cmd.replace("FASTA",FASTA).replace("universe",universe)
    os.system(carveme_cmd_mod)

    move_dmnd_cmd = "mv *.tsv ./dmnd_intermediates"
    os.system(move_dmnd_cmd)

def list_all_modules(Modules_directory):
    os.chdir(Modules_directory)
    list_all_mod = []
    all_mod = sorted(os.listdir())
    for module in all_mod:
        if module.endswith(".txt"):
            list_all_mod.append(module)
    return list_all_mod

def searchKeggShort(list_all_mod, Modules_directory, knumber):
    ''' Returns the list of Module_ids (Mxxxxx) in which a given KO is found '''
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
    '''
    Return a dictionary with KO as keys and reactions (R) as values, as per Module file.
    ------------------------------------------
    INPUT: KEGG Module flat file
    OUTPUT: dictionary manyXmany(?) of KEGG KO and KEGG Rs
    '''

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
    '''
    From a tabular data file that connects ModelSEED (old?) rxns to other onthologies,
    link KEGG RN to BiGG IDs
    & KEGG RN to ModelSEED IDs
    ------------------------------------------
    INPUT: reaction_DB (tabular data)
    OUTPUT: 1) dictionary manyXmany of KEGG Rs and BiGG IDs;
            2) dictionary oneXmany of KEGG Rs and ModelSEED IDs
    '''

    #TODO: reaction_DB could be updated with new releases of SEED and the rest

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
    v_KeggR_unique = []
    v_KeggR_unique_in_DB = []
    for Rnumber in list_of_unique_Rnumbers:
        v_KeggR_unique.append(Rnumber)
        if Rnumber in DB_Kegg_reactions.keys():
            v_KeggR_unique_in_DB.append(Rnumber)

    return v_KeggR_unique_in_DB

def bigg_nonredundant(v_KeggR_unique_in_DB, DB_Kegg_reactions):
    v_bigg_nonredundant = []
    for Rnumber in v_KeggR_unique_in_DB:
        for subdict in DB_Kegg_reactions[Rnumber]:
            for v_bigg in subdict.values():
                for single_bigg in v_bigg:
                    if single_bigg not in v_bigg_nonredundant:
                        v_bigg_nonredundant.append(single_bigg)

    return v_bigg_nonredundant

def modelseed_nonredundant(v_KeggR_unique_in_DB, DB_Kegg_reactions):
    v_modelseed_nonredundant = []
    for Rnumber in v_KeggR_unique_in_DB:
        for single_rxn in DB_Kegg_reactions[Rnumber]:
            if single_rxn not in v_modelseed_nonredundant:
                v_modelseed_nonredundant.append(single_rxn)

    return v_modelseed_nonredundant

def bigg_gapfill_absent_in_model(bigg_nonredundant_file, model):
    '''
    
    '''
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
    cmd_line = "curl --silent 'http://bigg.ucsd.edu/api/v2/universal/reactions/"
    el = single_reac
    command = cmd_line+el+"'"
    file_ift = str(el)+".ift"
    os.system(command+" > "+file_ift)

    return file_ift

def api_file_reorder(file_ift, verbose = False):
    with open (file_ift) as f:
        for line in f.readlines():
            text = "".join(line)
        #text cleaning
            text2 = text.replace(",", ",\n").replace(":",":\t").replace("&#8652;","<->")
            file_txt = str(file_ift).replace(".ift", ".txt")
            #files_txt.append(file_txt)
            g = open(file_txt, "w")
            g.write(text2)
            g.close()
            remove = "rm "
            os.system(remove+file_ift)

            if verbose:
                print("Info from file "+file_ift+" reordered")
    return 1

def api_file_reorder2(file_ift, verbose = False):
    with open (file_ift) as f:
        for line in f.readlines():
            text = "".join(line)
        #text cleaning
            text2 = text.replace(",", ",\n").replace(":",":\t").replace("&#8652;","<->")
            file_txt = str(file_ift).replace(".ift", ".txt")
            g = open(file_txt, "w")
            g.write(text2)
            g.close()
            remove = "rm "
            os.system(remove+file_ift)

            if verbose:
                print("Info from file "+file_ift+" reordered")

    return file_txt

def get_string(file_txt, total_strings, verbose = False, remove_intermediate = True):
    #extract reaction string and format it in Reframed format
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
            line = line.strip().replace("{", "") #TODO: check last-add: .replace("{", "")
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
        remove_cmd = "rm "
        os.system(remove_cmd+file_txt)

def get_stringa_error(file_txt, error_strings, verbose = False):
    error_strings.append(file_txt)
    if verbose:
        print(str(file_txt)+" string NOT PASSED") #debug

def retry_genes(file_txt, verbose = False):
    cmd_line = "curl --silent 'http://bigg.ucsd.edu/api/v2/search?query=GENE&search_type=genes'" #TODO: remove auto-quiet
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

            remove="rm "
            os.system(remove+file_ift2)
            return
        else:
            cmd_line2 = "curl --silent 'http://bigg.ucsd.edu/api/v2/models/MODEL/genes/BIGG'" #TODO: remove auto-quiet
            command3 = cmd_line2.replace("MODEL", model_bigg_id).replace("BIGG", bigg_id)
            file_ift3 = str(el)+".ift3"
            os.system(command3+" > "+file_ift3)

    remove="rm "
    os.system(remove+file_ift2)

    with open(file_ift3) as f:
        v = f.readline().replace(",","\n")
        new_single_reac = v.split('reactions": [{"bigg_id": "')[1].split('"\n')[0]
    new_ift = curl_bigg_reaction(new_single_reac)
    new_txt = api_file_reorder2(new_ift, verbose = verbose)
    try:
        get_string(new_txt, total_strings, verbose = verbose)
        error_strings.remove(file_txt)
        remove="rm "
        os.system(remove+file_ift3)
    except:
        return

def metab_change_names(model_new, verbose = False):
    # find the exact metabolite name and replace the one added via get_string()
    comm_list = []
    for metab in model_new.metabolites.values():
        metab_mod = str(metab).replace("M_", "", 1).replace("_c", "").replace("_e", "").replace("_p", "").replace("_m", "").replace("_n", "").replace("_x", "").replace("_h", "").replace("_g", "")
        if metab_mod in old_new_names.keys():
            command = "model_new.metabolites."+str(metab)+".name "+"="+' "'+old_new_names[metab_mod]+'"'
            comm_list.append(command)
            exec(command)

    if verbose:
        print(str(len(comm_list))+" commands used for METABOLITES name updates") #debug

    return comm_list

def check_name_changes(model_new):
    noncambiati=[]
    for x in model_new.metabolites.values():
        if str(x).startswith("M_"):
            noncambiati.append(x)

    if len(noncambiati) == 0:
        result = "All METABOLITES names were updated"
    if len(noncambiati) != 0:
        result = "ATTENTION Not all METABOLITES names were updated"

    print(result)

def reac_change_names(model_new, verbose = False):
    # find the exact reaction name and replace the one added via get_string()
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

    return comm_list

def old_new_names_dict(file):
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

def fixed_modules_of_interest(dir_base, module_file): # POSSIBILITY: CHOOSE MODULES FOR KO-RN CONNECTION (default: the ones in hmm.py)
    '''
    Point out MODULES from a fixed list as in the command instructions (.instruction file)
    '''
    MODofinterest = []

    os.chdir(dir_base)
    with open(module_file) as f:
        for line in f.readlines():
            MOD = line.strip()
            MODofinterest.append(MOD)
    return MODofinterest

def onbm_modules_of_interest(fasta_id, oneBM_modules_dir): # POSSIBILITY: CHOOSE MODULES FOR KO-RN CONNECTION (default: the ones in hmm.py)
    '''
    Point out MODULES with 1 block missing as in the command instructions (.instruction file)
    '''
    MODofinterest = []

    os.chdir(oneBM_modules_dir)
    if fasta_id+"_"+module_file in os.listdir():
        with open(fasta_id+"_"+module_file) as f:
            for line in f.readlines():
                MOD = line.strip()
                MODofinterest.append(MOD)
    return MODofinterest

def KOs_with_HMM_hits(hmm_hits_dir, fastakohits):
    '''
    Point out KOs from modules identified with hmm.py HITS. 
    ''' # formerly: those missing 1 block
    KOhits = []

    os.chdir(hmm_hits_dir)
    with open(fastakohits) as f:
        for line in f.readlines()[1:]:
            KO = line.strip().split("\t")[0]
            KOhits.append(KO)
    return KOhits

def Modules_KOhits_connection(KOhits, MODofinterest):
    '''
    Connect Modules with KOs missing, via dictionary.
    '''
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
    '''
    Connect Modules' names with the Modules flatfiles, via dictionary.
    '''
    modulesXflat = {}

    for module_txt in list_all_mod:
        for module in MODofinterestXKOhits.keys():
            if module_txt.startswith(module):
                modulesXflat.update({module:module_txt})

    return modulesXflat

def total_R_from_KOhits(MODofinterestXKOhits, modulesXflat, Modules_directory):
    '''
    Point out RN from KOs, scanning in Modules flatfiles.
    '''
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
    '''
    Save a log file of the output - FORMAT: 1 BiGG per line.
    '''
    os.chdir(gapfill_report_directory)
    f = open("bigg_log_"+fasta_id+".txt", "w")
    for single_bigg in bigg_nonredundant:
        f.write(single_bigg+"\n")
    f.close()

    print("COMPLETE BiGG logging "+fasta_id) #debug

def log_modelseed_nr(modelseed_nonredundant, fasta_id, gapfill_report_directory):
    '''
    Save a log file of the output - FORMAT: 1 SEED per line.
    '''
    os.chdir(gapfill_report_directory)
    f = open("seed_log_"+fasta_id+".txt", "w")
    for single_rxn in modelseed_nonredundant:
        f.write(single_rxn+"\n")
    f.close()

    print("COMPLETE ModelSEED logging "+fasta_id) #debug

def reframed_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = False):
    '''
    
    '''
    datetoday = str(datetime.datetime.now())[:10]
    # MODEL IO VIA REFRAMED
    os.chdir(model_directory)
    for file in sorted(os.listdir()):
        if file.endswith(".xml") and file.startswith(fasta_id):
            os.chdir(model_directory)
            model = load_cbmodel(file)

    # '''Task: saves reaction to add (from bigg_nonredundant() ) in a list'''
            os.chdir(gapfill_report_directory)
            bigg_gapfill = bigg_gapfill_absent_in_model("bigg_log_"+fasta_id+".txt", model)

    # '''Task: create a directory to store downloaded reaction from BiGG API'''
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
                with Pool(processes=6) as p: # POSSIBILITY: change number of processes/threads (now = 6)
                    p.map(curl_bigg_reaction, bigg_gapfill) #TO TAKE REACTIONS FROM SAVED LIST
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

    # '''Task: adds the reactions to a copy of the model, without those having any non-prokaryotic compartments'''
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

    # '''Task:  saves a LOGfile with the reactions ACTUALLY ADDED to the model'''
            os.chdir(gapfill_report_directory)
            f = open(""+model.id+"_added_reactions.txt", "w")
            for reac in added_reacs:
                f.write(reac+"\n")
            f.close()

    # '''Task:  checks namespace quality for metabolites and saves the new gap-filled model'''
            if verbose:
                metab_change_names(model_new, verbose = verbose)
                reac_change_names(model_new, verbose = verbose)
                check_name_changes(model_new)

            os.chdir(gapfilled_model_directory)
            save_cbmodel(model_new, model_new.id+"_KEGGadd_"+datetoday+".xml", flavor = 'bigg')

def cobrapy_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = False): #TODO: implement
    '''
    
    '''
    pass

def recap_addition(fasta_id, gapfill_report_directory, old_new_names_R):
    '''
    Stores RECAP informations for actual reaction addition in EACH model.
    It gives a final tabular output of the entire process.
    ------------------------------------------------------
    INPUT:  MAG/Genome of interest name - its contigs FASTA file should be in "genomes_directory" / input_directory
    OUTPUT: Table of added reactions strings.
    '''
    #TODO: enable single MAG/Genome output
    import datetime
    datetoday = str(datetime.datetime.now())[:10]

    os.chdir(gapfill_report_directory)
    f = open("recap_gapfill_"+datetoday+".tsv", "a")
    header = "MAG\tSTRING\tNAME\n"
    if f.tell() == 0:
        f.write(header)

    for file in sorted(os.listdir()):
        if file.endswith("_added_reactions.txt") and fasta_id in file:
            with open(file) as g:
                MAG = file[:-20]                                             #TODO: it works, this way, but better to avoid slicing
                for line in g.readlines():
                    STRING = line.strip()
                    reac_id = STRING.split(": ")[0].replace("R_", "")
                    NAME = old_new_names_R[reac_id]
                    if not STRING == "":
                        f.write(MAG+"\t"+STRING+"\t"+NAME+"\n")
    f.close()

if __name__ == "__main__":

    dir_base = os.getcwd()+"/"
    instruction_file = "genomes.instruction"
    module_file = "module_file.instruction"

    ####################
    # directories   DB #
    ####################
    Modules_directory = dir_base+"KEGG_MODULES_ref/"
    kkfiles_directory = Modules_directory+"file_test_kk/"

    ####################
    # directories  KMC #
    ####################
    KAnnotation_directory = dir_base+"KEGG_mappings/"

    ####################
    # directories  HMM #
    ####################
    hmm_hits_dir = dir_base+"HMM_HITS/"
    dir_genomes = dir_base+"genomes/"
    oneBM_modules_dir = dir_base+"oneBM_modules/"

    ####################
    # directories GSMM #
    ####################

    gapfill_report_directory = dir_base+"report_gapfill/"
    bigg_api = dir_base+"biggapi_download/"
    DB_directory = dir_base+"DB/"
    model_directory = dir_base+"models/"
    gapfilled_model_directory = dir_base+"models_gapfilled/"
    de_novo_model_directory = dir_base+"de_novo_models/"

########################################################################################

    parser = argparse.ArgumentParser(
        description="Genome-scale model reaction-addition after KEGG Modules Completeness evaluation and HMM-evidence for given genomes.")

    parser.add_argument('MAG_genome_FASTA',
                        help='''Run GSMM operation on the genome with the indicated fasta ID (no FASTA extension).''')
    parser.add_argument('-r','--reaction_addition_method', required=True,
                        choices=["reframed","cobrapy"],
                        help='''Selection of model I/O program for reaction addition.''')
    parser.add_argument('-d','--ontology_database_selection', required=True,
                        choices=["bigg","modelseed"],
                        help='''Selection of reaction ontology for gap-fill addition in GSMM.''') #TODO: implement

    parser.add_argument('--de_novo_gsmm', action ="store_true",
                        help='''Use all KEGG Modules missing 1 block for the HMM-based check (only KOs missing in the indicated MAG/Genome).''')

    parser.add_argument('--onebm_modules_list', action ="store_true",
                        help='''Use all KEGG Modules missing 1 block for the HMM-based check (only KOs missing in the indicated MAG/Genome).''')
    parser.add_argument('--fixed_modules_list', action ="store_true",
                        help='''Use a fixed list of KEGG Modules to use for the HMM-based check (only KOs missing in the indicated MAG/Genome).''') #TODO: check implementation
    parser.add_argument('-I','--path_input',
                        help='''Absolute path to input file(s) FOLDER.''', default = KAnnotation_directory)
    parser.add_argument('-O','--path_output',
                        help='''Absolute path to ouput file(s) FOLDER.''', default = dir_base)
    parser.add_argument('--log', action ="store_true",
                        help='''Store reaction-addition commands / info as plain text files.''') #TODO: implement
    parser.add_argument('-v','--verbose', action ="store_true",
                        help='''Print more informations - for debug or log.''')
    args = parser.parse_args()

########################################################################################

    #### SET DATA INFO
    list_all_mod = list_all_modules(Modules_directory)

    os.chdir(DB_directory)
    if args.ontology_database_selection == "bigg":
        DB_KEGG_RN = KEGG_BiGG_SEED_RN_dict("reactions_DB.tsv", DB_directory, ontology="BiGG")
        old_new_names = old_new_names_dict("metabolites_names_from_id_bigg.tsv")
        old_new_names_R = old_new_names_reac_dict("reactions_names_from_id_bigg.tsv")

    if args.ontology_database_selection == "modelseed":                                  #TODO: implement
        DB_KEGG_RN = KEGG_BiGG_SEED_RN_dict("reactions_DB.tsv", DB_directory, ontology="SEED")#TODO: implement
        old_new_names = old_new_names_dict("metabolites_names_from_id_seed.tsv")      #TODO: implement
        old_new_names_R = old_new_names_reac_dict("reactions_names_from_id_seed.tsv") #TODO: implement

########################################################################################

    ### READ INSTRUCTIONS
    os.chdir(dir_base)
    FASTA = str(args.MAG_genome_FASTA)
    with open(instruction_file) as f:
        for line in f.readlines()[1:]:
            if line.startswith(FASTA):
                line = line.strip().split("\t")
                for file in os.listdir(dir_genomes):
                    if file.startswith(FASTA):
                        fasta_genome = dir_genomes+file

                fasta_id = line[0]
                universe = line[2]
                klist_file = FASTA+".klist"
                fastakohits = FASTA+"_HMM_hits.txt"

    #### OPERATE SINGLE FUNCTIONS
    print("+++++START "+fasta_id)

    if args.de_novo_GSMM:
        build_de_novo_GSMM(FASTA, fasta_genome, de_novo_model_directory, hmm_hits_dir)

    else:
        if args.fixed_modules_list:
            MODofinterest = fixed_modules_of_interest(dir_base, module_file)
        if args.onebm_modules_list:
            MODofinterest = onbm_modules_of_interest(fasta_id, oneBM_modules_dir)
        KOhits = KOs_with_HMM_hits(hmm_hits_dir, fastakohits)
        MODofinterestXKOhits = Modules_KOhits_connection(KOhits, MODofinterest)
        modulesXflat = modules_flat_files_connection(list_all_mod, MODofinterestXKOhits)

        Rtotali_KOhits = total_R_from_KOhits(MODofinterestXKOhits, modulesXflat, Modules_directory)

        KEGG_R_to_add = keggR_in_DB(Rtotali_KOhits, DB_KEGG_RN)

# BiGG procedure
        if args.ontology_database_selection == "bigg":
            bigg_nonredundant = bigg_nonredundant(KEGG_R_to_add, DB_KEGG_RN)
            log_bigg_nr(bigg_nonredundant, fasta_id, gapfill_report_directory)

    # print out MAX ADDITION
            if args.verbose:
                os.chdir(gapfill_report_directory)
                for file in sorted(os.listdir()):
                    if file.startswith("bigg_log_"+fasta_id):
                        with open(file) as f:
                            print(fasta_id+" max reaction addition\t"+str(len(f.readlines())))

    #### REACTION ADDITION
            if args.reaction_addition_method == "reframed":
                if args.verbose:
                    reframed_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = True)
                else:
                    reframed_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = False)

            if args.reaction_addition_method == "cobrapy":
                if args.verbose:
                    cobrapy_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = True) #TODO
                else:
                    cobrapy_reaction_addition(fasta_id, model_directory, gapfill_report_directory, bigg_api, verbose = False) #TODO

# ModelSEED procedure
        elif args.ontology_database_selection == "modelseed":
            modelseed_nonredundant = modelseed_nonredundant(KEGG_R_to_add, DB_KEGG_RN)
            log_modelseed_nr(modelseed_nonredundant, fasta_id, gapfill_report_directory)

    # print out MAX ADDITION
            if args.verbose:
                os.chdir(gapfill_report_directory)
                for file in sorted(os.listdir()):
                    if file.startswith("seed_log_"+fasta_id):
                        with open(file) as f:
                            print(fasta_id+" max reaction addition\t"+str(len(f.readlines())))

    #### RECAP
    recap_addition(fasta_id, gapfill_report_directory, old_new_names_R)

print("END "+fasta_id+"\n")
