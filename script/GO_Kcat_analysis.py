import pandas as pd
import re
import os
import numpy as np
import json
from goatools import obo_parser
from pyprobar import bar, probar
from concurrent.futures import ThreadPoolExecutor
import plotly.graph_objects as go
import signal
import cobra

# 定义超时处理函数
def timeout_handler(signum, frame):
    raise TimeoutError("pFBA calculation timed out")

def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)
    
def calculate_depth(obo_path,go_term_id):
    go = obo_parser.GODag(obo_path)
    if go_term_id in go:
        go_term_info = go[go_term_id]
        
        return go_term_info.level

    else:
        print(f'GO Term with ID {go_term_id} is obsolete or not found.')
        
        return 0

def calculate_range_mean(value):
    if '-' in value:
        # 如果是范围，则计算范围的均值
        start, end = map(float, value.split('-'))
        return (start + end) / 2
    else:
        return float(value)

def calculate_mean(lst):
    kcat_values = list(map(float, lst))
    total = sum(kcat_values)
    count = len(kcat_values)
    mean = total / count
    return mean

def extract_kcat_GO_terms_from_protein_info(protein_info, GO_kcat_outfile):
    GO_pattern = r'GO:\d+'
    
    # 检查文件是否存在，如果不存在，则添加表头
    if not os.path.isfile(GO_kcat_outfile) or os.path.getsize(GO_kcat_outfile) == 0:
        with open(GO_kcat_outfile, 'w') as file:
            file.write("GO_term,Organism,Kcat\n")
            
    with open(GO_kcat_outfile, 'a') as file:
        for pkb_id, protein in protein_info.items():
            protein_kcat_values = []
            for substrate, substrate_kcat_info in protein['Properties']['Kcat'].items():
                for value in substrate_kcat_info:
                    protein_kcat_values.append(value['value'])
            if len(protein_kcat_values)>0:
                for entry, entry_info in protein['Consist'].items():
                    if 'Gene Ontology (GO)' in entry_info.keys():
                        Gene_Ontology = entry_info['Gene Ontology (GO)']
                        Gene_terms = re.findall(GO_pattern, Gene_Ontology)

                        for Gene_term in Gene_terms:
                            kcat_value = (';').join(protein_kcat_values)
                            org = '['+protein['Organism (ID)']+']'+protein['Organism']
                            file.write(Gene_term + ',' + str(org) + ',' + str(kcat_value) + '\n')
                            

def extract_kcat_GO_terms_from_protein_info_return_df(protein_info,GO_term_kcat_realtion_df,GO_term_kcat_realtion_df_index):
    GO_pattern = r'GO:\d+'
    for pkb_id,protein in protein_info.items():
        protein_kcat_values = []
        for substrate,substrate_kcat_info in protein['Properties']['Kcat'].items():
            for value in substrate_kcat_info:
                protein_kcat_values.append(value['value'])
                
        for entry,entry_info in protein['Consist'].items():
            if 'Gene Ontology (GO)' in entry_info.keys():
                Gene_Ontology = entry_info['Gene Ontology (GO)']
                Gene_terms = re.findall(GO_pattern, Gene_Ontology)
                for Gene_term in Gene_terms:
                    # print(Gene_term)
                    kcat_values = []
                    GO_term_kcat_realtion_df.loc[GO_term_kcat_realtion_df_index,'GO_term'] = Gene_term
                    GO_term_kcat_realtion_df.loc[GO_term_kcat_realtion_df_index,'Kcat'] = (';').join(protein_kcat_values)
                    
                    GO_term_kcat_realtion_df_index+=1
    
    return GO_term_kcat_realtion_df,GO_term_kcat_realtion_df_index

def get_GO_term_mean_kcat_from_GO_kcat_file(protein_GO_kcat_file,obo_path,go_term_id):
    GO_term_kcat_realtion_df = pd.read_csv(protein_GO_kcat_file)
    GO_term_kcat_realtion_df = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['Kcat'].notnull()] 
    
    go = obo_parser.GODag(obo_path)
    if go_term_id in go:
        go_term = go[go_term_id]
        
        # print("GO Term ID:", go_term_id)
        # print("GO Term Name:", go_term.name)
        # print("Namespace:", go_term.namespace)

        child_term_kcat_sum = 0

        # Get the child terms (direct children)
        #print("\nChild Terms:")
        child_term_len = len(go_term.get_all_children())
        if child_term_len:
            for child_term in go_term.get_all_children():
                #print(child_term)

                if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]) > 0:
                    kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]['Kcat'].values[0]
                    # print(kcat_values)
                    # kcat_list = [float(value) for value in kcat_values.split(';')]
                    # kcat_list = [float(value) for value in kcat_values.split(';') if '-' not in value]

                    kcat_list = [calculate_range_mean(value) for value in kcat_values.split(';')]

                    # print("kcat_list:", kcat_list)

                    if kcat_list:
                        kcat_mean = calculate_mean(kcat_list)

                        # print('kcat_mean:',kcat_mean)
                        child_term_kcat_sum += kcat_mean

            # print(f'{go_term_id}——MEAN kcat:',child_term_kcat_sum/child_term_len)

            return child_term_kcat_sum/child_term_len
        
        else:
            return '-'
        
    else:
        print(f'GO Term with ID {go_term_id} is obsolete or not found.')
        return '-'
    
def get_kcat_from_GO(protein_GO_kcat_file,obo_path,go_term_id):
    # print(go_term_id)
    GO_term_kcat_realtion_df = pd.read_csv(protein_GO_kcat_file)
    GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['Kcat'].notnull()] 
    go_term_kcat_sum = 0

    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]) > 0:
        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]['Kcat'].values[0]

        # kcat_list = [float(value) for value in kcat_values.split(';')]
        # kcat_list = [float(value) for value in kcat_values.split(';') if '-' not in value]
        kcat_list = [calculate_range_mean(value) for value in kcat_values.split(';')]
        
        kcat_mean = calculate_mean(kcat_list)

        print(f'{go_term_id}——MEAN kcat:',kcat_mean)
        
        return kcat_mean
    
    else:
        go_term_mean_kcat = get_GO_term_mean_kcat_from_GO_kcat_file(protein_GO_kcat_file,obo_path,go_term_id)
        print(f'{go_term_id}——go_term_mean_kcat',go_term_mean_kcat)
        
        return go_term_mean_kcat

def get_kcat_from_GO_by_orgname_old(org_name,GO_kcat_file,obo_path,go_term_id):
    org_name_file = org_name.replace(' ','_')
    GO_term_kcat_df = pd.read_csv(GO_kcat_file)
    # Use str.contains on the 'Organism' column
    org_GO_term = GO_term_kcat_df[GO_term_kcat_df['Organism'].str.contains(org_name)]
    org_GO_term[org_GO_term['Kcat'].notnull()] 
    org_GO_kcat_file = f'./analysis/{org_name_file}_GO_kcat.csv'
    org_GO_term.to_csv(org_GO_kcat_file, header=True, index=False, mode='w')

    go_term_mean_kcat = get_kcat_from_GO(org_GO_kcat_file,obo_path,go_term_id)
    
    return [go_term_mean_kcat,org_GO_kcat_file]

def get_kcat_from_GO_by_orgid_old(org_id,GO_kcat_file,obo_path,go_term_id):
    GO_term_kcat_df = pd.read_csv(GO_kcat_file)

    GO_term_kcat_df['Organism_ID'] = GO_term_kcat_df['Organism'].str.extract(r'\[(\d+)\]')
    org_GO_term = GO_term_kcat_df[GO_term_kcat_df['Organism_ID'] == org_id][['GO_term','Organism','Kcat']]
    org_GO_kcat_file = f'./analysis/{org_id}_GO_kcat.csv'
    org_GO_term.to_csv(org_GO_kcat_file, header=True, index=False, mode='w')

    go_term_mean_kcat = get_kcat_from_GO(org_GO_kcat_file,obo_path,go_term_id)
    return go_term_mean_kcat,org_GO_kcat_file

def get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,go_term_id):
    GO_term_kcat_df = pd.read_csv(GO_kcat_file)

    GO_term_kcat_df['Organism_ID'] = GO_term_kcat_df['Organism'].str.extract(r'\[(\d+)\]')
    GO_term_kcat_realtion_df = GO_term_kcat_df[GO_term_kcat_df['Organism_ID'] == org_id][['GO_term','Organism','Kcat']]
    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]) > 0:
        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]['Kcat'].values
        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]
        kcat_mean = calculate_mean(kcat_list)
        #print(f'{go_term_id}——MEAN kcat:',kcat_mean)
        return kcat_mean
    else:
        #go = obo_parser.GODag(obo_path)
        if go_term_id in go:
            go_term = go[go_term_id]
            
            # print("GO Term ID:", go_term_id)
            # print("GO Term Name:", go_term.name)
            # print("Namespace:", go_term.namespace)

            child_term_kcat_sum = 0

            # Get the child terms (direct children)
            #print("\nChild Terms:")
            child_term_len = len(go_term.get_all_children())
            if child_term_len:
                for child_term in go_term.get_all_children():
                    #print(child_term)
                    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]) > 0:
                        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]['Kcat'].values
                        # print(kcat_values)
                        # kcat_list = [float(value) for value in kcat_values.split(';')]
                        # kcat_list = [float(value) for value in kcat_values.split(';') if '-' not in value]

                        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]

                        # print("kcat_list:", kcat_list)

                        if kcat_list:
                            kcat_mean = calculate_mean(kcat_list)

                            # print('kcat_mean:',kcat_mean)
                            child_term_kcat_sum += kcat_mean

                # print(f'{go_term_id}——MEAN kcat:',child_term_kcat_sum/child_term_len)
                return child_term_kcat_sum/child_term_len
            else:
                return '-'
        else:
            print(f'GO Term with ID {go_term_id} is obsolete or not found.')
            return '-'

def get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,go_term_id):
    GO_term_kcat_df = pd.read_csv(GO_kcat_file)
    # Use str.contains on the 'Organism' column
    GO_term_kcat_realtion_df = GO_term_kcat_df[GO_term_kcat_df['Organism'].str.contains(org_name)]
    GO_term_kcat_realtion_df = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['Kcat'].notnull()] 
    #print(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id])
    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]) > 0:
        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]['Kcat'].values
        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]
        kcat_mean = calculate_mean(kcat_list)
        #print(f'{go_term_id}——MEAN kcat:',kcat_mean)
        return kcat_mean
    
    else:
        #go = obo_parser.GODag(obo_path)
        if go_term_id in go:
            go_term = go[go_term_id]
            
            # print("GO Term ID:", go_term_id)
            # print("GO Term Name:", go_term.name)
            # print("Namespace:", go_term.namespace)

            child_term_kcat_sum = 0

            # Get the child terms (direct children)
            #print("\nChild Terms:")
            child_term_len = len(go_term.get_all_children())
            if child_term_len:
                for child_term in go_term.get_all_children():
                    #print(child_term)
                    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]) > 0:
                        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]['Kcat'].values
                        # print(kcat_values)
                        # kcat_list = [float(value) for value in kcat_values.split(';')]
                        # kcat_list = [float(value) for value in kcat_values.split(';') if '-' not in value]

                        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]

                        # print("kcat_list:", kcat_list)

                        if kcat_list:
                            kcat_mean = calculate_mean(kcat_list)

                            # print('kcat_mean:',kcat_mean)
                            child_term_kcat_sum += kcat_mean

                # print(f'{go_term_id}——MEAN kcat:',child_term_kcat_sum/child_term_len)
                return child_term_kcat_sum/child_term_len
            else:
                return '-'
        else:
            print(f'GO Term with ID {go_term_id} is obsolete or not found.')
            return '-'

def get_go_term_mean_kcat_by_org(norm_model,ns2assc_org,org_type,org_name,org_id,GO_kcat_file,go):
    go_term_mean_kcat_dict = {} 
    # reaction-gene-ncbigene-GO-kcat
    for rea in probar(norm_model.reactions):  
        #print(rea.id)
        go_term_mean_kcat_dict[rea.id] = {}
        go_term_mean_kcat_dict[rea.id]['BP'] = {}
        go_term_mean_kcat_dict[rea.id]['CC'] = {}
        go_term_mean_kcat_dict[rea.id]['MF'] = {}    
        go_term_mean_kcat_dict[rea.id]['Total'] = {}  

        gene_list = rea._gene_reaction_rule.replace('(','').replace(')','').replace('and','').replace('or','').split() 
        for each_g in gene_list:
            try:
                norm_model.genes.get_by_id(each_g).annotation['ncbigene']
            except:
                print('Gene has no annatation!')
            else:
                g2ncbi = int(norm_model.genes.get_by_id(each_g).annotation['ncbigene'])
                
                #BP
                try:
                    ns2assc_org['BP'][g2ncbi]
                except:
                    go_term_mean_kcat_dict[rea.id]['BP'][each_g] = []
                else:
                    BP_GO_set = ns2assc_org['BP'][g2ncbi]
                    mean_kcat_list =[]
                    for eachgo in BP_GO_set:
                        if org_type =='org_name':
                            go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,eachgo)
                        elif org_type =='org_id':
                            go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,eachgo)
                        mean_kcat_list.append(go_term_mean_kcat)
                        #print(rea.id,each_g,eachgo,go_term_mean_kcat)
                    go_term_mean_kcat_dict[rea.id]['BP'][each_g] = mean_kcat_list
                #MF
                try:
                    ns2assc_org['MF'][g2ncbi]
                except:
                    go_term_mean_kcat_dict[rea.id]['MF'][each_g] = []
                else:
                    MF_GO_set = ns2assc_org['MF'][g2ncbi]
                    mean_kcat_list =[]
                    for eachgo in MF_GO_set:
                        if org_type =='org_name':
                            go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,eachgo)
                        elif org_type =='org_id':
                            go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,eachgo)
                        mean_kcat_list.append(go_term_mean_kcat)
                        #print(rea.id,each_g,eachgo,go_term_mean_kcat)
                    go_term_mean_kcat_dict[rea.id]['MF'][each_g] = mean_kcat_list
                #CC
                try:
                    ns2assc_org['CC'][g2ncbi]
                except:
                    go_term_mean_kcat_dict[rea.id]['CC'][each_g] = []
                else:
                    CC_GO_set = ns2assc_org['CC'][g2ncbi]
                    mean_kcat_list =[]
                    for eachgo in CC_GO_set:
                        if org_type =='org_name':
                            go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,eachgo)
                        elif org_type =='org_id':
                            go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,eachgo)
                        mean_kcat_list.append(go_term_mean_kcat)
                        #print(rea.id,each_g,eachgo,go_term_mean_kcat)
                    go_term_mean_kcat_dict[rea.id]['CC'][each_g] = mean_kcat_list
                go_term_mean_kcat_dict[rea.id]['Total'][each_g] = go_term_mean_kcat_dict[rea.id]['BP'][each_g]+go_term_mean_kcat_dict[rea.id]['MF'][each_g]+go_term_mean_kcat_dict[rea.id]['CC'][each_g]
        #break
    return go_term_mean_kcat_dict

#运行一直没结果。。。
def get_go_term_mean_kcat_by_org_concurrent(go_term_mean_kcat_dict,eachr,norm_model,ns2assc_org,org_type,org_name,org_id,GO_kcat_file,go):
    print(eachr)
    go_term_mean_kcat_dict[eachr] = {}
    go_term_mean_kcat_dict[eachr]['BP'] = {}
    go_term_mean_kcat_dict[eachr]['CC'] = {}
    go_term_mean_kcat_dict[eachr]['MF'] = {}    
    go_term_mean_kcat_dict[eachr]['Total'] = {}  
    rea = norm_model.reactions.get_by_id(eachr)
    gene_list = rea._gene_reaction_rule.replace('(','').replace(')','').replace('and','').replace('or','').split() 
    for each_g in gene_list:
        try:
            norm_model.genes.get_by_id(each_g).annotation['ncbigene']
        except:
            print('Gene has no annatation!')
        else:
            g2ncbi = int(norm_model.genes.get_by_id(each_g).annotation['ncbigene'])
            
            #BP
            try:
                ns2assc_org['BP'][g2ncbi]
            except:
                go_term_mean_kcat_dict[rea.id]['BP'][each_g] = []
            else:
                BP_GO_set = ns2assc_org['BP'][g2ncbi]
                mean_kcat_list =[]
                for eachgo in BP_GO_set:
                    if org_type =='org_name':
                        go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,eachgo)
                    elif org_type =='org_id':
                        go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,eachgo)
                    mean_kcat_list.append(go_term_mean_kcat)
                    #print(rea.id,each_g,eachgo,go_term_mean_kcat)
                go_term_mean_kcat_dict[rea.id]['BP'][each_g] = mean_kcat_list
            #MF
            try:
                ns2assc_org['MF'][g2ncbi]
            except:
                go_term_mean_kcat_dict[rea.id]['MF'][each_g] = []
            else:
                MF_GO_set = ns2assc_org['MF'][g2ncbi]
                mean_kcat_list =[]
                for eachgo in MF_GO_set:
                    if org_type =='org_name':
                        go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,eachgo)
                    elif org_type =='org_id':
                        go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,eachgo)
                    mean_kcat_list.append(go_term_mean_kcat)
                    #print(rea.id,each_g,eachgo,go_term_mean_kcat)
                go_term_mean_kcat_dict[rea.id]['MF'][each_g] = mean_kcat_list
            #CC
            try:
                ns2assc_org['CC'][g2ncbi]
            except:
                go_term_mean_kcat_dict[rea.id]['CC'][each_g] = []
            else:
                CC_GO_set = ns2assc_org['CC'][g2ncbi]
                mean_kcat_list =[]
                for eachgo in CC_GO_set:
                    if org_type =='org_name':
                        go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name,GO_kcat_file,go,eachgo)
                    elif org_type =='org_id':
                        go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id,GO_kcat_file,go,eachgo)
                    mean_kcat_list.append(go_term_mean_kcat)
                    #print(rea.id,each_g,eachgo,go_term_mean_kcat)
                go_term_mean_kcat_dict[rea.id]['CC'][each_g] = mean_kcat_list
            go_term_mean_kcat_dict[rea.id]['Total'][each_g] = go_term_mean_kcat_dict[rea.id]['BP'][each_g]+go_term_mean_kcat_dict[rea.id]['MF'][each_g]+go_term_mean_kcat_dict[rea.id]['CC'][each_g]

    return go_term_mean_kcat_dict

def get_go_term_mean_kcat_by_org_v2old(norm_model, ns2assc_org, org_type, org_name, org_id, GO_kcat_file, go):
    go_term_mean_kcat_dict = {}

    def get_mean_kcat_list(gene, go_type):
        mean_kcat_list = []
        try:
            ncbigene = int(norm_model.genes.get_by_id(gene).annotation['ncbigene'])
            go_set = ns2assc_org[go_type].get(ncbigene, set())
            for eachgo in go_set:
                if org_type == 'org_name':
                    go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name, GO_kcat_file, go, eachgo)
                elif org_type == 'org_id':
                    go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id, GO_kcat_file, go, eachgo)
                if go_term_mean_kcat !='-':
                    mean_kcat_list.append(go_term_mean_kcat)
        except:
            print('Error occurred while processing gene:', gene)
        return mean_kcat_list

    for rea in probar(norm_model.reactions):
        gene_list = rea._gene_reaction_rule.replace('(','').replace(')','').replace('and','').replace('or','').split() 
        go_term_mean_kcat_dict[rea.id] = {
            'BP': {gene: get_mean_kcat_list(gene, 'BP') for gene in gene_list },
            'MF': {gene: get_mean_kcat_list(gene, 'MF') for gene in gene_list },
            'CC': {gene: get_mean_kcat_list(gene, 'CC') for gene in gene_list }
        }
        go_term_mean_kcat_dict[rea.id]['Total'] = {
            gene: go_term_mean_kcat_dict[rea.id]['BP'][gene] +
                  go_term_mean_kcat_dict[rea.id]['MF'][gene] +
                  go_term_mean_kcat_dict[rea.id]['CC'][gene]
            for gene in go_term_mean_kcat_dict[rea.id]['BP']
        }

    return go_term_mean_kcat_dict

def get_go_term_mean_kcat_by_org_v2(norm_model, ns2assc_org, org_type, org_name, org_id, GO_kcat_file, go):
    go_term_mean_kcat_dict = {}

    def get_mean_kcat_list(gene, go_type):
        mean_kcat_list = []
        try:
            ncbigene = int(norm_model.genes.get_by_id(gene).annotation['ncbigene'])
            go_set = ns2assc_org[go_type].get(ncbigene, set())
            for eachgo in go_set:
                if org_type == 'org_name':
                    go_term_mean_kcat = get_kcat_from_GO_by_orgname(org_name, GO_kcat_file, go, eachgo)
                elif org_type == 'org_id':
                    go_term_mean_kcat = get_kcat_from_GO_by_orgid(org_id, GO_kcat_file, go, eachgo)
                if go_term_mean_kcat !='-':
                    mean_kcat_list.append(go_term_mean_kcat)
        except:
            print('Error occurred while processing gene:', gene)
        return mean_kcat_list

    def process_reaction(rea):
        gene_list = rea._gene_reaction_rule.replace('(','').replace(')','').replace('and','').replace('or','').split() 
        go_term_mean_kcat_dict[rea.id] = {
            'BP': {gene: get_mean_kcat_list(gene, 'BP') for gene in gene_list },
            'MF': {gene: get_mean_kcat_list(gene, 'MF') for gene in gene_list },
            'CC': {gene: get_mean_kcat_list(gene, 'CC') for gene in gene_list }
        }
        go_term_mean_kcat_dict[rea.id]['Total'] = {
            gene: go_term_mean_kcat_dict[rea.id]['BP'][gene] +
                  go_term_mean_kcat_dict[rea.id]['MF'][gene] +
                  go_term_mean_kcat_dict[rea.id]['CC'][gene]
            for gene in go_term_mean_kcat_dict[rea.id]['BP']
        }

    with ThreadPoolExecutor() as executor:
        executor.map(process_reaction, norm_model.reactions)

    return go_term_mean_kcat_dict

def process_data(go_term_mean_kcat_dict, go_type, clc_type):
    reaction_kcat_dict = {}
    for reaction, data in go_term_mean_kcat_dict.items():
        gene_values_clean = {gene: [value for value in values if value != '-'] for gene, values in data[go_type].items()}
        gene_values = {gene: values for gene, values in gene_values_clean.items() if values}
        if gene_values:
            values = [value for sublist in gene_values.values() for value in sublist]
            if clc_type == 'max':
                reaction_kcat_dict[reaction] = np.max(values)
            elif clc_type == 'mean':
                reaction_kcat_dict[reaction] = np.mean(values)
            elif clc_type == 'median':
                reaction_kcat_dict[reaction] = np.median(values)
    return reaction_kcat_dict

def get_kcat_from_GO(org_id,org_name,use_type,GO_kcat_file,go,go_term_id):
    GO_term_kcat_df = pd.read_csv(GO_kcat_file)

    if use_type == 'org_name':
        # Use str.contains on the 'Organism' column
        GO_term_kcat_realtion_df = GO_term_kcat_df[GO_term_kcat_df['Organism'].str.contains(org_name)]
        GO_term_kcat_realtion_df = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['Kcat'].notnull()] 
    elif use_type == 'org_id':
        GO_term_kcat_df['Organism_ID'] = GO_term_kcat_df['Organism'].str.extract(r'\[(\d+)\]')
        GO_term_kcat_realtion_df = GO_term_kcat_df[GO_term_kcat_df['Organism_ID'] == org_id][['GO_term','Organism','Kcat']]
        GO_term_kcat_realtion_df = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['Kcat'].notnull()] 
    else:
        GO_term_kcat_realtion_df = GO_term_kcat_df[GO_term_kcat_df['Kcat'].notnull()] 

    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]) > 0:
        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == go_term_id]['Kcat'].values
        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]
        kcat_mean = calculate_mean(kcat_list)
        #print(f'{go_term_id}——MEAN kcat:',kcat_mean)
        return kcat_mean
    else:
        #go = obo_parser.GODag(obo_path)
        if go_term_id in go:
            go_term = go[go_term_id]
            
            # print("GO Term ID:", go_term_id)
            # print("GO Term Name:", go_term.name)
            # print("Namespace:", go_term.namespace)

            child_term_kcat_sum = 0

            # Get the child terms (direct children)
            #print("\nChild Terms:")
            child_term_len = len(go_term.get_all_children())
            if child_term_len:
                for child_term in go_term.get_all_children():
                    #print(child_term)
                    if len(GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]) > 0:
                        kcat_values = GO_term_kcat_realtion_df[GO_term_kcat_realtion_df['GO_term'] == child_term]['Kcat'].values
                        # print(kcat_values)
                        # kcat_list = [float(value) for value in kcat_values.split(';')]
                        # kcat_list = [float(value) for value in kcat_values.split(';') if '-' not in value]

                        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]

                        # print("kcat_list:", kcat_list)

                        if kcat_list:
                            kcat_mean = calculate_mean(kcat_list)

                            # print('kcat_mean:',kcat_mean)
                            child_term_kcat_sum += kcat_mean

                # print(f'{go_term_id}——MEAN kcat:',child_term_kcat_sum/child_term_len)
                return child_term_kcat_sum/child_term_len
            else:
                return '-'
        else:
            print(f'GO Term with ID {go_term_id} is obsolete or not found.')
            return '-'

def get_go_term_mean_kcat_by_org_v3(norm_model, ns2assc_org, org_type, org_name, org_id, GO_kcat_file, go):
    go_term_mean_kcat_dict = {}

    def get_mean_kcat_list(gene, go_type):
        mean_kcat_list = []
        try:
            ncbigene = int(norm_model.genes.get_by_id(gene).annotation['ncbigene'])
            go_set = ns2assc_org[go_type].get(ncbigene, set())
            for eachgo in go_set:
                go_term_mean_kcat = get_kcat_from_GO(org_id,org_name,org_type,GO_kcat_file,go,eachgo)
                if go_term_mean_kcat !='-':
                    mean_kcat_list.append(go_term_mean_kcat)
        except:
            print('Error occurred while processing gene:', gene)
        return mean_kcat_list

    def process_reaction(rea):
        gene_list = rea._gene_reaction_rule.replace('(','').replace(')','').replace('and','').replace('or','').split() 
        go_term_mean_kcat_dict[rea.id] = {
            'BP': {gene: get_mean_kcat_list(gene, 'BP') for gene in gene_list },
            'MF': {gene: get_mean_kcat_list(gene, 'MF') for gene in gene_list },
            'CC': {gene: get_mean_kcat_list(gene, 'CC') for gene in gene_list }
        }
        go_term_mean_kcat_dict[rea.id]['Total'] = {
            gene: go_term_mean_kcat_dict[rea.id]['BP'][gene] +
                  go_term_mean_kcat_dict[rea.id]['MF'][gene] +
                  go_term_mean_kcat_dict[rea.id]['CC'][gene]
            for gene in go_term_mean_kcat_dict[rea.id]['BP']
        }

    with ThreadPoolExecutor() as executor:
        executor.map(process_reaction, norm_model.reactions)

    return go_term_mean_kcat_dict

def compare_json(file1, file2):
    data1 = load_json(file1)
    data2 = load_json(file2)
    
    return compare_dicts_sort(data1, data2)

def filter_json(data):
    if isinstance(data, dict):
        return {k: filter_json(v) for k, v in data.items() if v != "-"}
    elif isinstance(data, list):
        return [filter_json(v) for v in data if v != "-"]
    else:
        return data
    
def compare_dicts_sort(dict1, dict2):
    diff = {}
    
    # Check keys in dict1
    for key in dict1:
        if key not in dict2:
            diff[key] = ("Key only in dict1", dict1[key], None)
        elif dict1[key] != dict2[key]:
            if isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
                diff[key] = compare_dicts(dict1[key], dict2[key])
            else:
                # Sort lists before comparison
                if isinstance(dict1[key], list) and isinstance(dict2[key], list):
                    dict1[key].sort()
                    dict2[key].sort()
                diff[key] = ("Value differs", dict1[key], dict2[key])
    
    # Check keys in dict2
    for key in dict2:
        if key not in dict1:
            diff[key] = ("Key only in dict2", None, dict2[key])
    
    return diff

def compare_dicts_sort_formal(dict1, dict2):
    diff = {}
    
    # Check keys in dict1
    for key in dict1:
        if key not in dict2:
            diff[key] = ("Key only in dict1", round(dict1[key], 3), None)
        elif round(dict1[key], 3) != round(dict2[key], 3):
            if isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
                diff[key] = compare_dicts(dict1[key], dict2[key])
            else:
                # Sort lists before comparison
                if isinstance(dict1[key], list) and isinstance(dict2[key], list):
                    dict1[key].sort()
                    dict2[key].sort()
                diff[key] = ("Value differs", round(dict1[key], 3), round(dict2[key], 3))
    
    # Check keys in dict2
    for key in dict2:
        if key not in dict1:
            diff[key] = ("Key only in dict2", None, round(dict2[key], 3))
    
    return diff

def draw_boxplot_cbfig_diff_cond(growth_pfba_enzmodel,columns_to_calculate,diff_model_diff_substate_result_figfile):
    # Calculate errors
    errors = {
        column: abs(growth_pfba_enzmodel['EXP'] - growth_pfba_enzmodel[column]) / growth_pfba_enzmodel['EXP']
        for column in columns_to_calculate
    }

    # Create traces
    traces = [go.Box(y=errors[column], name=f'<b>{column}<b>') for column in errors]

    # Define layout
    layout = go.Layout(
        xaxis=dict(title=dict(font=dict(size=20, family='Times New Roman')), tickfont=dict(color='black', size=15, family='Times New Roman')),
        yaxis=dict(title=dict(text="<b>Estimation error<b>", font=dict(size=20, family='Times New Roman')), tickfont=dict(color='black', size=15, family='Times New Roman')),
        showlegend=False, width=1200, height=800
    )

    # Create figure
    fig = go.Figure(data=traces, layout=layout)
    fig.update_layout(xaxis_tickangle=45)

    # Write image and show figure
    fig.write_image(diff_model_diff_substate_result_figfile, scale=1, width=1200, height=800)
    fig.show()

def set_model_parameter_by_substrate(growth_model,substrate,ori_substrate_id_list,concentration):
    for eachsubid in ori_substrate_id_list:
        if re.search('_reverse',eachsubid):
            r_id_new=eachsubid.split('_reverse')[0]
            growth_model.reactions.get_by_id(eachsubid).bounds = (0, 0) 
            growth_model.reactions.get_by_id(r_id_new).bounds = (0, 0)  
        else:
            r_id_new=eachsubid+'_reverse'
            growth_model.reactions.get_by_id(eachsubid).bounds = (0, 0) 
            growth_model.reactions.get_by_id(r_id_new).bounds = (0, 0) 
            
    growth_model.reactions.get_by_id(substrate).bounds = (-concentration, 0)
    try:
        growth_model.reactions.get_by_id(substrate+'_reverse')
    except:
        pass
    else:
        growth_model.reactions.get_by_id(substrate+'_reverse').bounds = (0, 0)
        
    return growth_model

def draw_cbfig_diff_cond(growth_rate_diff_substrate_file,method):
    ECMpy_diff_substate_result_figfile ="./analysis/%s_diff_substate_result.png"%method
    growth_pfba_enzmodel = pd.read_csv(growth_rate_diff_substrate_file,index_col=0)
    #color_dim=np.random.randint(0,255,len(growth_pfba_enzmodel['EXP']))
    color_dim=['#eaff56','#fa8c35','#bce672','#ff0097','#ff7500','#ffc64b','#00bc12','#c0ebd7',\
            '#eacd76','#eedeb0','#21a675','#808080','#88ada6','#725e82','#3d3b4f','#a78e44',\
            '#3de1ad','#75664d','#ff2d51','#ffb3a7','#44cef6','#4b5cc4','#8d4bbb','#3b2e7e']
    
    #Create traces
    data1=[]
    i=0
    for eachmet in list(growth_pfba_enzmodel.index):  
        trace0 = go.Scatter(name=eachmet,x=[growth_pfba_enzmodel.loc[eachmet,'EXP']], y=[growth_pfba_enzmodel.loc[eachmet,method]], 
                            marker={'color': color_dim[i], 'size': 20}, mode='markers',showlegend=True)
        #marker={'color': color_dim[i], 'size': 20, 'symbol':i}
        data1.append(trace0)
        i=i+1

    trace1 = go.Scatter(x = [0,2],y = [0,2],mode = 'lines',marker = dict(color = 'rgb(127, 127, 127)'),showlegend=False)
    data1.append(trace1)

    layout = go.Layout(
        xaxis= dict(title=dict(text="<b>Measured growth rate (h<sup>-1</sup>)<b>", font=dict(size=20,family='Times New Roman')),range=[0,1],tickfont=dict(color= 'black',size=15,family='Times New Roman')),
        yaxis=dict(title=dict(text="<b>Predicted growth rate (h<sup>-1</sup>)<b>", font=dict(size=20,family='Times New Roman')),range=[0,2],tickfont=dict(color= 'black',size=15,family='Times New Roman')),
        showlegend=True,width=1200,height=800)
    fig = go.Figure(data=data1, layout=layout)
    fig.update_layout(
        legend=dict(
            title_font_family="Times New Roman",  # 图例标题字体
            font=dict(  # 图例字体
                family="Times New Roman",
                size=16,
                color="black"  # 颜色：红色
            )
        )
    )
    fig.write_image(ECMpy_diff_substate_result_figfile, scale=1, width=1200, height=800)
    fig.show()

def calculate_normalized_error(data, exp_fluxes_error):
    fluxes_error = (np.sum(list(map(lambda num:num*num, data)))/(len(data)-1))** 0.5
    return fluxes_error / exp_fluxes_error

def print_statistics(name, data, normalized_error):
    arr_mean = np.mean(data)
    arr_var = np.var(data)
    arr_std = np.std(data, ddof=1)
    print("%s mean value is: %.6f" % (name, arr_mean))
    print("%s var is: %.6f" % (name, arr_var))
    print("%s sd is: %.6f" % (name, arr_std))
    print("%s normalized error is: %.6f\n" % (name, normalized_error))

def calculate_model_substrate_pfba_old(model,substrate,concentration,obj,result_dict):
    with model as growth_model: 
        result = {}
        # change original substrate in model
        [ori_obj_id,ori_substrate_id_list,ori_sub_concentration,ori_ATPM] = get_model_substrate_obj(growth_model)
        growth_model = set_model_parameter_by_substrate(growth_model, substrate, ori_substrate_id_list, concentration)
        
        # 设置超时时间（秒）
        timeout_seconds = 10
        
        # 注册信号处理程序
        signal.signal(signal.SIGALRM, timeout_handler)
        
        # 设置计算超时时间
        signal.alarm(timeout_seconds)
        
        try:
            cobra.flux_analysis.pfba(growth_model)
            # 如果成功完成，取消超时
            signal.alarm(0)
        except TimeoutError as e:
            print(e)
            result_dict[substrate] = 0
            # 取消超时
            signal.alarm(0)
        else:
            pfba_solution = cobra.flux_analysis.pfba(growth_model)
            result_dict[substrate] = pfba_solution.fluxes[obj]
    return result_dict

# 放到GO_Kcat_analysis里面调用，result_dict未空
def calculate_model_substrate_pfba(model, substrate, concentration, obj, result_dict):
    with model as growth_model: 
        result = {}
        # change original substrate in model
        [ori_obj_id, ori_substrate_id_list, ori_sub_concentration, ori_ATPM] = get_model_substrate_obj(growth_model)
        growth_model = set_model_parameter_by_substrate(growth_model, substrate, ori_substrate_id_list, concentration)
        
        while True:  # 无限循环，直到满足条件才跳出
            start_time = time.time()  # 每次循环开始时重新设置计时器
            try:
                pfba_solution = cobra.flux_analysis.pfba(growth_model)
                result_dict[substrate] = pfba_solution.fluxes[obj]
                break  # 成功计算完成后跳出循环
            except TimeoutError as e:
                print(e)
                result_dict[substrate] = 0
                break  # 超时后立即跳出循环
            if time.time() - start_time > 10:  # 判断计算时间是否超过10秒
                break  # 超时后立即跳出循环