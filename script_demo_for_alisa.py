import os
import pandas as pd
import requests, sys
import pandas as pd
import json
from pandas import json_normalize
import time
import requests, sys
import json
import numpy as np
import sys
import argparse
import pprint
import optparse,os,sys,subprocess,time


'''
parser = argparse.ArgumentParser(description='A test program.')

parser.add_argument("-i", "--input_file", help="Reads in chr, pos, ref, alt file.")
parser.add_argument("-o", "--output_file", help="Creates file in this new location.")

args = parser.parse_args()

input_file= str(args.input_file)
output_file= str(args.output_file)
'''
    
def setup_df(chrom_col_name, pos_col_name, rsid_col_name, ref_col_name, alt_col_name, input_file, col_sep):
    if input_file[-4:] ==".vcf":

        with open(input_file, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        test= pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
        chrom_col_name="CHROM"
        pos_col_name="POS"
        ref_col_name="REF"
        alt_col_name="ALT"
        rsid_col_name="ID"
        full_df=test

    else:
        test= pd.read_csv(input_file, sep= col_sep)


    full_df= test

    cols_to_insert=["indel?","VEP alleles (ref/alt1/alt2/...)", "VEP ref", "VEP alt", "VEP worst consequence",
                    "Impact","Symbol", "Returned consequence types","Unique biotypes",
                    "All consequence terms", "Exon", "Intron","FathmmXF noncoding prediction", "FathmmXF coding prediction"]
    
    if rsid_col_name== "": #if there's no rsid col, add it to the col to insert
        cols_to_insert.insert(0,"rsid")
        rsid_col_name="rsid"


    counter=0
    insertion_index= full_df.shape[1]
    for col_name in cols_to_insert:
        if col_name not in full_df.columns:
            full_df.insert(loc=insertion_index +counter, column=col_name, value="")
        counter +=1
    index_col_header=0

    return full_df
    #print(str(full_df.loc[full_df.index[3], "RSID"]))
    #making a package "import __"

    
    
def parse_and_insert_to_df(chrom_col_name, pos_col_name, rsid_col_name, ref_col_name, alt_col_name, constructor, full_df): #manaully check on https://rest.ensembl.org/vep/human/id/rs1225841234?appris=1&CADD=1
    global gene_names_all
    global possible_keys
    global count
    server = "https://rest.ensembl.org"
    #extra_annotations= "AncestralAllele",
    ext = "/vep/human/region?AncestralAllele=1&transcript_version=1&numbers=1&vcf_string=1"+plugins_for_url
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    formatted= '{ "variants" : [%s]}' % (constructor)
    #print(formatted)
    tic= time.perf_counter()
    r = requests.post(server+ext, headers=headers, data=formatted)
    toc=time.perf_counter()

    print(f"Request for {count-1} in {toc - tic:0.4f} seconds")
    tic= time.perf_counter()

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    query_return=eval(repr(decoded))
    #print(query_return)
    search_key=set()

    for index_of_dict_obj, dictionary_obj in enumerate(query_return):

        gene_names=set()
        

        unique_vep_types= set()
        gene_sym_worst_conseq= set()

        impact= set()
        symbol_gene= set()
        ENSG_id= set()
        ENST_id= set()
        biotype_worst_consequence= set()
        unique_motif_conseq= set()
        all_consequence_terms= set()
        exon= set()
        intron= set()
        motif_info=[]


        #print(dictionary_obj)
        #print('-'*50)
        input_query= dictionary_obj['input']
        input_query=input_query.split()
        #print(input_query)
        chrom_query= input_query[0]
        pos_query= input_query[1]
        rsid_query= input_query[2]
        ref_query= input_query[3]
        alt_query= input_query[4]
        multiallelic= False
        for key in dictionary_obj.keys(): #to detect if unexpected fields are found
            possible_keys.add(key)

        #vep_alleles= dictionary_obj["allele_string"]


        indices=full_df.loc[(full_df[chrom_col_name].astype(str) == chrom_query)& 
                        (full_df[pos_col_name].astype(str) == pos_query) & 
                        ((full_df[rsid_col_name].astype(str) == "") | (full_df[rsid_col_name].astype(str) == rsid_query)) &
                        (full_df[ref_col_name].astype(str) == ref_query) & 
                        (full_df[alt_col_name].astype(str) == alt_query)]

        indexes_returned=np.array(indices.index.tolist())

        index_in_df= indexes_returned[0] #index of the first row since we are looking at unique inputs

        worst_consequences_list= dictionary_obj["most_severe_consequence"]
        returned_consequence_types = set(key for key in dictionary_obj.keys() if key.endswith("_consequences"))


        #queries through the possible terms
        consequence_fields= {'regulatory_feature_consequences' : "REGULATORY", 'transcript_consequences' : "TRANSCRIPT",
                             'motif_feature_consequences' : "MOTIF", 'intergenic_consequences' : "INTERGENIC"
                            }
        #relevant_motif_fields= ['high_inf_pos', 'motif_name', 'motif_pos', 'motif_score_change']
        for conseq, value in consequence_fields.items(): 
            biotype_by_conseq= set()
            ENSx_by_conseq= set()
            ENSx_key=[]
            if conseq in dictionary_obj:
                for term in dictionary_obj[conseq]:
                    temp_motif=[]
                    
                    available_fields=term.keys()
                    ENSx_key= [i for i in available_fields if conseq[:-12] in i] #searches for the ENSx key in term
                    for consequence_term in term["consequence_terms"]:
                        all_consequence_terms.add(consequence_term)
                    if 'biotype' in term:
                        unique_vep_types.add(f"{value}: {term['biotype']}")
                        biotype_by_conseq.add(term['biotype'])
                   
                        
                        
                    if len(ENSx_key) >=1 and ENSx_key[0] in term:
                        ENSx_by_conseq.add(term[ENSx_key[0]])
                        if conseq == 'motif_feature_consequences':
                            temp_motif.append(term[ENSx_key[0]])
                            if 'motif_score_change' in term:
                                temp_motif.append(term['motif_score_change'])
                            if 'high_inf_pos' in term:
                                temp_motif.append(term['high_inf_pos'])

                    if len(temp_motif) >=1:    
                        motif_info.append(temp_motif)
                            
                        
                        
                    ancestral_allele=""
                    if worst_consequences_list in term["consequence_terms"]:
                        if 'biotype' in term:
                            biotype_worst_consequence.add(term['biotype'])
                        if 'gene_id' in term:
                            ENSG_id.add(term['gene_id'])
                        if 'transcript_id' in term:
                            ENST_id.add(term['transcript_id'])
                        if 'impact' in term:
                            impact.add(term['impact'])
                        if 'exon' in  term:
                            exon.add(term['exon'])
                        if 'intron' in term:
                            intron.add(term['intron'])
                        if 'variant_allele' in term:
                            variant_allele= term['variant_allele']
                        if 'aa' in term:
                            ancestral_allele= term['aa']

                        

                        if 'gene_symbol' in term:
                            gene_sym_worst_conseq.add(term['gene_symbol'])
                            
                            symbol_gene.add(term['gene_symbol'])
                            gene_names_all.add(term['gene_symbol'])
                        for plugin_name, plugin_value in numeric_cols.items():

                            if plugin_value[0] is not int and plugin_value[0] in available_fields: #no plugins needed
                                full_df.at[index_in_df, plugin_value[0]] = term[plugin_value[0]]

                            elif len(plugin_value) == 1 and any(plugin_name.lower() in sub for sub in available_fields) and plugin_name not in search_key: #keys aren't known
                                print(plugin_name)
                                search_key.add(plugin_name)
                                        
                                print("-"*50) 
                                pprint.pprint(term)
                            if (len(plugin_value) > 1 and plugin_value[1] in available_fields) or plugin_name == "dbNSFP":
                                if plugin_name== "mutfunc": #FINISHED I THINK
                                    mutfunc_mode= term[plugin_value[1]]
                                    for mutfunc_type, mutfunc_dict in mutfunc_mode.items():
                                        for mutfunc_col, mutfunc_value in mutfunc_dict.items():
                                            full_df.at[index_in_df, f"mutfunc_{mutfunc_type}_{mutfunc_col}"] = mutfunc_value
                                            
                                elif plugin_name == "GeneSplicer":
                                    genesplicer_multiple_sites=term[plugin_value[1]].split(",")
                                    genesplicer_parse={"state":0, # {col_name : index after "\" split}
                                                       "type": 1,
                                                       "coordinates":2,
                                                       "confidence":3,
                                                       "score":4,
                                                      }
                                    for genesplicer_count, site in enumerate(genesplicer_multiple_sites):
                                        for genesplicer_key, genesplicer_value in genesplicer_parse.items():
                                            full_df.at[index_in_df, f"genesplicer_site{genesplicer_count}_{genesplicer_key}"] = site.split("/")[genesplicer_value]

                                elif plugin_name == "LoF":
                                    for lof_key in plugin_value[1:]:
                                        if lof_key == "lof_info":
                                            lof_info_str= term[lof_key].split(",")
                                            for lof_each_info in lof_info_str:
                                                lof_col_and_val = lof_each_info.split(":") #Example PERCENTILE:0.765765765765766
                                                full_df.at[index_in_df, f"lof_{lof_col_and_val[0]}"] = lof_col_and_val[1]
                                        else:
                                            if lof_key in available_fields:
                                                full_df.at[index_in_df, lof_key] = term[lof_key]
                                elif plugin_name == "DisGeNET":
                                    full_df.at[index_in_df, 'disgenet_dict']= str(term[plugin_value[1]])


                                elif plugin_name == "dbNSFP":
                                    dbNSFP_fields= plugin_value[0].lower().split(',')
                                    for dbNSFP_annotation in dbNSFP_fields:
                                        if dbNSFP_annotation in term:
                                            full_df.at[index_in_df, f"dbNSFP_{dbNSFP_annotation}"] = term[dbNSFP_annotation]

                                            
                                elif plugin_name == "SpliceAI":   
                                    for spliceai_key, spliceai_value in term[plugin_value[1]].items():
                                        full_df.at[index_in_df, f"spliceai_{spliceai_key}"] = spliceai_value
                                else:
                                    for normal_parse in plugin_value[1:]: #parses through the keys and inserts into df if they exist
                                        if normal_parse in available_fields:
                                            full_df.at[index_in_df, normal_parse] = term[normal_parse]
                #print(index_in_df)
                if conseq == 'motif_feature_consequences':
                    full_df.loc[[index_in_df], 'motif_info'] = pd.Series([motif_info], index=full_df.index[[index_in_df]])
                
                elif ENSx_key:
                    #print(f"{ENSx_key[0]} -> {ENSx_by_conseq}")
                    full_df.loc[[index_in_df], ENSx_key[0]] = pd.Series([list(ENSx_by_conseq)], index=full_df.index[[index_in_df]])
                full_df.loc[[index_in_df], f"{value}_biotypes"] = pd.Series([list(biotype_by_conseq)], index=full_df.index[[index_in_df]])





        #formatting the variables
        #filters out duplicate gene_symbols picked up


        rsid_from_colocated_var= ""
        #if there's no colocated variants OR if colocated variants doesn't have rsid
        if 'colocated_variants' not in dictionary_obj or not dictionary_obj['colocated_variants'][-1]['id'].startswith("rs"):
            #print("Colocated_variants not in")
            #pprint.pprint(dictionary_obj)
            #print('__'*50)

            vep_alleles= dictionary_obj["allele_string"]

        else:
            rsid_info= dictionary_obj['colocated_variants'][-1]
            vep_alleles= rsid_info['allele_string']
            rsid_from_colocated_var= rsid_info['id']

            if full_df.loc[index_in_df, rsid_col_name] == "":
                full_df.at[index_in_df, rsid_col_name] = rsid_from_colocated_var

            if 'frequencies' in rsid_info:
                for var_allele in rsid_info['frequencies'].keys():
                    for pop_study, freq in rsid_info['frequencies'][var_allele].items():
                        full_df.at[index_in_df, pop_study] = freq
        
        
        vep_allele_list=vep_alleles.split("/")



        #str_unique_types= list(unique_vep_types)                               


        #for allele in allele_freq_dict:
        symbol_gene=', '.join(symbol_gene)
        impact=', '.join(impact)

        #adding data to df
        column_name_val ={'VEP worst consequence': worst_consequences_list, 'Impact': impact, 'Symbol': symbol_gene,
                          'worst conseq gene id': list(ENSG_id), 'worst conseq transcript id': list(ENST_id), 
                          'Returned consequence types': list(returned_consequence_types), 'Exon': list(exon), 'Intron': list(intron),
                          'VEP alleles (ref/alt1/alt2/...)': vep_alleles, 'VEP ref': ancestral_allele, 'VEP alt': variant_allele,
                          'Unique biotypes': list(unique_vep_types), 'All consequence terms': list(all_consequence_terms)
                         }

        for col_name, value in column_name_val.items():
            try:
                if isinstance(value, float) or len(value) != 0: #if it's not empty
                    full_df.at[index_in_df, col_name] = value
            except Exception:
                full_df.loc[[index_in_df], col_name] = pd.Series([col_name], index=full_df.index[[index_in_df]])
    toc=time.perf_counter()
    print(f"Parsing the request took a total of {toc - tic:0.4f} seconds")
    return full_df




count=1
numeric_cols={#"Blosum62":[1, "blosum62"], #plugin_name: [bool_val, key_in_query] #this plugin doesn't work see 2-16-23 notes
              #"UTRAnnotator":[1, "Existing_uORFs", "Existing_OutOfFrame_oORFs", "Existing_InFrame_oORFs"], #not found yet
              
              "CADD":[1,"cadd_phred","cadd_raw"],
              "EVE":[1, "eve_score", "eve_class"],
              "IntAct":[1], #key found, but may not contain numeric info
              "LoF":[1, "lof", "lof_filter", "lof_flag", "lof_info"],
              #"Phenotypes":[1], #this plugin looks to be very resource heavy
              "GeneSplicer":[1, "genesplicer"],
              "dbscSNV":[1, "rf_score", "ada_score"], #key found
             
              #"dbNSFP":["LRT_pred,MutationTaster_pred,REVEL_rankscore,PrimateAI_score", 1], #What predictions should we retreive?
              "DisGeNET":[1, "disgenet"],
              "SpliceAI":[2, "spliceai"],
              "MaxEntScan":[1, "maxentscan_alt", "maxentscan_diff", "maxentscan_ref"], 
              #"mutfunc":[1, "mutfunc"],
              "polyphen":["polyphen_score"],
              "sift": ["sift_score"]
             }


plugins_for_url=""
for plugin_name, value in numeric_cols.items():
    if type(value[0]) is int or plugin_name== "dbNSFP":
        plugins_for_url+=f"{plugin_name}={value[0]}&"

plugins_for_url=plugins_for_url[:-1]
gene_names_all=set()
possible_keys=set()
expected_keys= ['colocated_variants', 'assembly_name', 'end', 'allele_string', 'start', 
                'seq_region_name', 'vcf_string', 'input', 'regulatory_feature_consequences', 
                'motif_feature_consequences', 'most_severe_consequence', 
                'strand', 'id', 'transcript_consequences', 'intergenic_consequences']


def create_batch_query(chrom_col_name, pos_col_name, rsid_col_name, ref_col_name, alt_col_name, full_df, output_file):
    global count
    duplicate_tracker=set()
    constructor=""
    original_full_df= full_df
    for index_for_last, row in original_full_df.iterrows():#loops through each entry

        chrom=str(row[chrom_col_name])
        target_pos= str(row[pos_col_name]) #initializes parameters
        rsid=str(row[rsid_col_name])
        if type(rsid) is float or rsid =="": #some files have null and some don't have rsid col -> defaults to ""
            rsid= "nan"
        ref_allele= str(row[ref_col_name])
        alt_allele= str(row[alt_col_name])

        full_df.at[index_for_last,'indel?'] = False
        if len(ref_allele) > 1 or len(alt_allele) >1: #detects indels and changes status, but still queries using API
            full_df.at[index_for_last,'indel?'] = True

        str_for_constructor = f'"{chrom} {target_pos} {rsid} {ref_allele} {alt_allele} . . .", '


        if f"{target_pos} {ref_allele} {alt_allele}" not in duplicate_tracker: #POSITION DUPLICATE FILTER
            duplicate_tracker.add(f"{target_pos} {ref_allele} {alt_allele}")
            constructor += str_for_constructor
            #print(constructor)
            #print()
            count+=1

        if count % 200 == 0 or index_for_last == original_full_df.shape[0]-1:
            constructor = constructor[0:-2]
            #print(len(duplicate_tracker))
            #print(row)
            full_df= parse_and_insert_to_df(chrom_col_name, pos_col_name, rsid_col_name, ref_col_name, alt_col_name, constructor, full_df)
            #print('-'*50)

            constructor=""
            count=1

    if not (set(possible_keys).issubset(set(expected_keys))):
        print(set(possible_keys) ^ set(expected_keys))            

    full_df.to_csv(str(output_file)+ ".csv", index=False)
    
    
    return full_df, gene_names_all


def write_geneHancer_search(full_df, chrom_col_name, pos_col_name, output_file):
    with open(f"{output_file}GH_instructions.txt", 'w') as f:
        f.write("For GeneHancer: copy and paste as 'defined regon' into https://genome.ucsc.edu/cgi-bin/hgTables\n")
        f.write("select group=Regulation and track=GeneHancer and output=tsv\n")
        for chrom_num in full_df[chrom_col_name].unique():
            subsetted_by_chrom= full_df[full_df[chrom_col_name] == chrom_num]
            if chrom_num == 23:
                chrom_num="X"
            
            
            f.write(f"chr{chrom_num} {subsetted_by_chrom[pos_col_name].min()} {subsetted_by_chrom[pos_col_name].max()}\n")
            
def parse_genehancer(chrom_col_name, pos_col_name, full_df, genehancer_file, output_file):

    GeneHancer_df=pd.read_csv(genehancer_file, sep='\t')

    #Gene to GHID, verified using chromosome and position
    for index, row in full_df.iterrows():#loops through each entry
        filter1= GeneHancer_df[GeneHancer_df["#chrom"] == f"chr{row[chrom_col_name]}"]
        filter2= filter1[filter1["chromStart"] <= int(row[pos_col_name])]
        filter3= filter2[filter2["chromEnd"] >= int(row[pos_col_name])]
        if filter3.empty:
            GH_id_only= "GHID not found"
        else:
            GH_id_only=filter3['name'].drop_duplicates().item()
            full_df.at[index,'GeneHancer Score'] = filter3['score'].item()
        full_df.at[index,'GeneHancer ID'] = GH_id_only

    full_df.to_csv(str(output_file)+ ".csv", index= False)

def write_fathmm_search(chrom_col_name, pos_col_name, ref_col_name, alt_col_name, full_df, output_file):
    with open(f"{output_file}_FATHMM.txt", 'w') as f:
        f.write("for FATHMM-XF: copy and paste into https://fathmm.biocompute.org.uk/fathmm-xf/ -> Download and manually save (ctrl + s) as tsv file for input file in part2\n") 
        set_to_insert= set()
        for index, row in full_df.iterrows():#loops through each entry to create a query for copy/paste
            if len(row[ref_col_name]) == 1 and len(row[alt_col_name]) == 1:
                chrom=str(row[chrom_col_name])
                position= str(row[pos_col_name])
                ref_allele= str(row[ref_col_name])
                alt_allele= str(row[alt_col_name])
                set_to_insert.add(f'{chrom},{position},{ref_allele},{alt_allele}')
        for string_to_insert in set_to_insert:
            f.write(string_to_insert + '\n')
                #https://fathmm.biocompute.org.uk/fathmm-xf/cgi-bin/results.cgi?session=8a9b9fe7-12c1-4b81-9111-2e03d224f444
def parse_fathmm(chrom_col_name, pos_col_name, rsid_col_name, ref_col_name, alt_col_name, full_df, fathmm_file, output_file):
    test= pd.read_csv(fathmm_file, sep="\t", names=['Chromosome','Position','Ref. Base','Mutant Base','Coding Score','Non-Coding Score','Warning'])
    test=test.iloc[2: , :]

    XF_df= test
    already_queried_xf=[]
    for index, row in full_df.iterrows():#loops through each entry
        chrom=row[chrom_col_name]
        position= row[pos_col_name]
        ref_allele= row[ref_col_name]
        alt_allele= row[alt_col_name]
        if len(ref_allele) == 1 and len(alt_allele) == 1 and position not in already_queried_xf:
            already_queried_xf.append(position)
            XF_row= XF_df[(XF_df['Chromosome'] == str(chrom)) &
                          (XF_df['Position'] == str(position)) &
                          (XF_df['Ref. Base']== ref_allele) &
                          (XF_df['Mutant Base']== alt_allele)]

            non_coding_val=np.array(XF_row['Non-Coding Score'].tolist())
            coding_val=np.array(XF_row['Coding Score'].tolist())

            if len(non_coding_val) > 0 and len(coding_val) > 0:
                if non_coding_val[0] != "--":
                    full_df.at[index,'FathmmXF noncoding prediction'] = non_coding_val[0]
                elif coding_val[0] != "--":
                    full_df.at[index,'FathmmXF coding prediction'] = coding_val[0]

            else:
                XF_val= ""
                full_df.at[index,'FathmmXF noncoding prediction'] = XF_val
                full_df.at[index,'FathmmXF coding prediction'] = XF_val

    full_df.to_csv(str(output_file)+ ".csv", index=False)
    return full_df
def parse_ABC_131_biosample(chrom_col_name, pos_col_name, ABC_input, full_df, output_file):
    ABC_hg38_output = pd.read_csv(ABC_input)
    for index, row in full_df.iterrows():
        chrom= row[chrom_col_name]
        pos= int(row[pos_col_name])
        ABC_row= ABC_hg38_output[(ABC_hg38_output["hg38_chr"] == f"chr{chrom}") &
                                 (ABC_hg38_output["hg38_start"] < pos) &
                                 (ABC_hg38_output["hg38_end"] > pos)]
        if len(ABC_row) > 0:
            full_df.at[index, "ABC Unique TargetGene"] = ABC_row.TargetGene.unique()
            full_df.at[index, "ABC Unique CellTypes"] = ABC_row.CellType.unique()
            
            full_df.at[index, "TargetGene with highest ABC score"] = f'[{ABC_row.loc[ABC_row["ABC.Score"].idxmax()]["TargetGene"]}, {ABC_row["ABC.Score"].max()}]'
    full_df.to_csv(str(output_file)+ ".csv", index=False)
        

def parse_111_bing_ren(chrom_col_name, pos_col_name, ABC_input, full_df, output_file):
    db=ABC_input
    for index, row in full_df.iterrows():
        (chrom,pos) = f"chr{row[chrom_col_name]}", int(row[pos_col_name])
        unique_genes= set()
        unique_cells=set()
        cmd         = "tabix %s %s:%d-%d" % (db, chrom, pos-1, pos+1)
        print(cmd)
        proc        = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, err = proc.communicate()

        result=result.decode()
        print(result)
        highest_ABC_score=0
        highest_score_geneName=""
        if len(result) > 0:
            for record in result.split("\n"):
                if not record:
                    continue
                parts  = record.strip().split("\t")
                if float(parts[-4]) > float(highest_ABC_score):
                    highest_ABC_score=parts[-4]
                    highest_score_geneName=parts[-3].split(":")[0]
                    full_df.at[index, "ABC Bing Ren Highest Score TargetGene"] = f"[{highest_score_geneName}, {highest_ABC_score}]"

                
                unique_genes.add(parts[-3].split(":")[0])
                unique_cells.add(parts[-1])



        if len(unique_genes) > 0:
            full_df.at[index, "ABC Bing Ren Unique TargetGene"] = str(unique_genes)
            full_df.at[index, "ABC Bing Ren Unique cells"] = str(unique_cells)
    full_df.to_csv(str(output_file)+ ".csv", index=False)

def write_vcf(chrom_col_name, pos_col_name, rsid_col_name, ref_col_name, alt_col_name, full_df, output_file):
    with open(f"{output_file}.vcf", 'w') as f: #creates pseudo vcf file

        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\n")

        for index, row in full_df.iterrows():
            chrom=str(row[chrom_col_name])
            position= str(row[pos_col_name])
            rsid= str(row[rsid_col_name])
            ref_allele= str(row[ref_col_name])
            alt_allele= str(row[alt_col_name])
            f.write(f'{chrom}\t{position}\t{rsid}\t{ref_allele}\t{alt_allele}\n') #no ID column
