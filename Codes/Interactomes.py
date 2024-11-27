import pandas as pd
from unipressed import IdMappingClient
import time
import copy
import pickle
from modules import check_condition, check_condition_reverse
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import matplotlib
import requests

#Load saved variables
with open('.\Variables\ensp_to_uniprot.pkl', 'rb') as f:
    ensp_to_uniprot = pickle.load(f)
with open(r'.\Variables\gname_to_pfam.pkl', 'rb') as f:
    gname_to_pfam = pickle.load(f)
with open('.\Variables\gname_to_uniprot.pkl', 'rb') as f:
    gname_to_uniprot = pickle.load(f)
with open(r'.\Variables\uniprot_to_pfam.pkl', 'rb') as f:
    uniprot_to_pfam = pickle.load(f)
with open(r'.\Variables\non_kinase_genes.pkl', 'rb') as f:
    non_kinase_genes = pickle.load(f)
net_prop = pd.read_csv(r'..\Results\network_properties.tsv', sep="\t", index_col=0)

##Read the overlap results from "Parsing_the input_files.py"
overlaps = pd.read_csv('..\Results\overlaps_phosphorylation.tsv', sep='\t')
overlaps['KinaseID'] = overlaps['KinaseID'].apply(lambda x: str(x[2:-3])) #Format
overlaps_NetworKIN = pd.read_csv(r'..\Results\NetworKINvsPhospho.ELM.tsv', sep='\t')
overlaps_NetworKIN['KinaseID'] = overlaps_NetworKIN['KinaseID'].apply(lambda x: str(x[2:-3])) #Format
overlaps_GPS = pd.read_csv('..\Results\GPSvsPhosphoSitePlus.tsv', sep='\t')

##Interactomes databases
#STRING
STRING = pd.read_csv('..\Data\9606.protein.links.v11.5.txt', sep=" ")
# region kinase-kinase interactions
# kinase_kinase = pd.read_csv('..\Data\string_interactions.tsv', sep="\t")
# #Filter out only kinase-substrate interactions
# #Remove kinase-kinase interactions
# STRING_subset = STRING[~STRING[['protein1', 'protein2']].apply(lambda x: tuple(x), axis=1).isin(kinase_kinase[['node1_string_id', 'node2_string_id']].apply(lambda x: tuple(x), axis=1))]
# STRING_subset = STRING_subset.copy() #To overcome "SettingWithCopyWarning" message
# endregion
#Remove the prefix "9606."
STRING_subset = STRING.copy() #Added
STRING_subset['protein1'] = STRING_subset['protein1'].str.replace('9606.', '', regex=True)
STRING_subset['protein2'] = STRING_subset['protein2'].str.replace('9606.', '', regex=True)

#Convert Ensembl protein IDs to Uniprot IDs
# region ensp_to_uniprot
# ensp = list(set(STRING_subset['protein1'].tolist()+STRING_subset['protein2'].tolist()))
# request = IdMappingClient.submit(
#     source="Ensembl_Protein", dest="UniProtKB", ids=ensp
# )
# time.sleep(15.0)
# converted = list(request.each_result())
# #There are 18777 unique Ensembl protein IDs, out of which only 18364 were converted to UniProt IDs.
# ensp = list(map(lambda x: x['from'], converted))
# uniprot = list(map(lambda x: x['to'], converted))
# for i in range(len(ensp)):
#     ensp_to_uniprot[ensp[i]]= uniprot[i]
# with open('.\Variables\ensp_to_uniprot.pkl', 'wb') as f:
#     pickle.dump(ensp_to_uniprot, f)
# endregion
unipressed = copy.deepcopy(STRING_subset)
unipressed[['protein1', 'protein2']] = unipressed[['protein1', 'protein2']].applymap(lambda x: ensp_to_uniprot[x] if x in ensp_to_uniprot else pd.NA)
unipressed.dropna(inplace=True)
unipressed.drop_duplicates(inplace=True)

# region BioMart
# #'ensp_to_uniprot' contains IDs from both BioMart and unipressed, so this comparison is obsolete.
# #BioMart
# biomart = pd.read_csv('..\Data\mart_export.txt', sep='\t')
# biomart.drop_duplicates(inplace=True)
# biomart.dropna(inplace=True)
# # biomart['Protein stable ID'].isin(ensp).sum()
# ensp_biomart = biomart['Protein stable ID'][biomart['Protein stable ID'].isin(ensp)].tolist()
#
# #BioMart vs Unipressed
# ensp_unipressed = list(ensp_to_uniprot.keys())
# diff = list(set(ensp_biomart) - set(ensp_unipressed))
# # diff = list(set(ensp_unipressed) - set(ensp_biomart))
#
# biomart_dict = dict(zip(biomart['Protein stable ID'], biomart['UniProtKB/Swiss-Prot ID']))
# BioMart = copy.deepcopy(STRING_format)
# BioMart[['protein1', 'protein2']] = BioMart[['protein1', 'protein2']].applymap(lambda x: biomart_dict[x] if x in biomart_dict else pd.NA)
# BioMart.dropna(inplace=True)
#
# merge1 = pd.merge(overlaps, BioMart, left_on=['SubstrateID', 'KinaseID'], right_on=['protein1', 'protein2'])
# merge2 = pd.merge(overlaps, BioMart, left_on=['SubstrateID', 'KinaseID'], right_on=['protein2', 'protein1'])
# overlaps_biomart = pd.concat([merge1, merge2])
# overlaps_biomart = overlaps_biomart.drop(['protein1', 'protein2'], axis=1)
# overlaps_biomart.drop_duplicates(inplace=True)
#
# merge1 = pd.merge(overlaps, unipressed, left_on=['SubstrateID', 'KinaseID'], right_on=['protein1', 'protein2'])
# merge2 = pd.merge(overlaps, unipressed, left_on=['SubstrateID', 'KinaseID'], right_on=['protein2', 'protein1'])
# overlaps_unipressed = pd.concat([merge1, merge2])
# overlaps_unipressed = overlaps_unipressed.drop(['protein1', 'protein2'], axis=1)
# overlaps_unipressed.drop_duplicates(inplace=True)
# endregion

#Map uniprot IDs to protein families to filter out those which have "kinase"
# region uniprot_to_pfam
# uniprot= list(set(unipressed['protein1'].tolist()+unipressed['protein2'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# for id in uniprot[17000:]:
#         query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
#         response = requests.get(query_url)
#         protein_name = response.text.split('\n')[1]
#         if "kinase" in protein_name.lower():
#             uniprot_to_pfam[id] = protein_name
# with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_pfam, f)
# endregion
#Retain phosphorylation interactions
kinases = uniprot_to_pfam.keys()
STRING_A = unipressed[unipressed['protein1'].isin(kinases)]
STRING_B = unipressed[unipressed['protein2'].isin(kinases)]
STRING_subset = pd.concat([STRING_A, STRING_B])
del unipressed
# region kinase-kinase interactions
# #Remove kinase-kinase interactions
# STRING_format = STRING_subset[~((STRING_subset['protein1'].isin(kinases)) & (STRING_subset['protein2'].isin(kinases)))]
# endregion
#No duplicates
#STRING contains bidirectional interactions. Drop them to get exact kinase-substrate interactions
STRING_format = STRING_subset.copy() #Added
STRING_format = STRING_format[~STRING_format.apply(lambda x: tuple(sorted(x[['protein1', 'protein2']])), axis=1).duplicated(keep='first')]
# Check for autophosphorylation
# kk = STRING_format[STRING_format['protein1']==STRING_format['protein2']] # No interactions
# STRING_format.to_csv(r'.\Variables\STRING_format.tsv', sep="\t", index=False)

#Overlaps
#STRING vs prediction tools
merge1 = pd.merge(overlaps, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein1', 'protein2'])
merge2 = pd.merge(overlaps, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein2', 'protein1'])
STRING_overlaps = pd.concat([merge1, merge2])
STRING_overlaps = STRING_overlaps.drop(['protein1', 'protein2'], axis=1)
STRING_overlaps.to_csv(r'..\Results\STRINGvsPredictiontools.tsv', sep="\t", index=False)

#STRING vs NetworKIN and Phospho.ELM database
merge1 = pd.merge(overlaps_NetworKIN, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein1', 'protein2'])
merge2 = pd.merge(overlaps_NetworKIN, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein2', 'protein1'])
STRING_NetworKIN = pd.concat([merge1, merge2])
STRING_NetworKIN = STRING_NetworKIN.drop(['protein1', 'protein2'], axis=1)
STRING_NetworKIN.to_csv(r'..\Results\STRINGvsOverlapsNetworKIN.tsv', sep="\t", index=False)

#STRING vs GPS and PhosphoSitePlus database
merge1 = pd.merge(overlaps_GPS, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein1', 'protein2'])
merge2 = pd.merge(overlaps_GPS, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein2', 'protein1'])
STRING_GPS = pd.concat([merge1, merge2])
STRING_GPS = STRING_GPS.drop(['protein1', 'protein2'], axis=1)
STRING_GPS.to_csv(r'..\Results\STRINGvsOverlapsGPS.tsv', sep="\t", index=False)

#BioGRID
BioGRID = pd.read_csv(r'..\Data\BIOGRID-ORGANISM-Homo_sapiens-4.4.219.tab3.txt', sep="\t", dtype='unicode')
#dtype='unicode' to overcome DtypeWarning: Columns (1,2,18) have mixed types. Specify dtype option on import or set low_memory=False.
BioGRID_subset = BioGRID.loc[BioGRID['Organism Name Interactor A'].eq('Homo sapiens') & BioGRID['Organism Name Interactor B'].eq('Homo sapiens')]
BioGRID_subset = BioGRID_subset.iloc[:,[7,8,23,26]]
#There are missing UniProt IDs, but there are no missing gene names.
BioGRID_subset.columns = BioGRID_subset.columns[:2].tolist()+['protein1', 'protein2']
#Map gene names to protein families to filter out those which have "kinase"
# region gname_to_pfam
# gene_names= list(set(BioGRID_subset['Official Symbol Interactor A'].tolist()+BioGRID_subset['Official Symbol Interactor B'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# gene_names = [x for x in gene_names if x not in gname_to_pfam.keys()]
# for gene_name in gene_names[19000:]:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_name = response.text.split('\n')[1]
#     if "kinase" in protein_name.lower():
#         gname_to_pfam[gene_name] = protein_name
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
# non_kinase_genes = [x for x in gene_names if x not in gname_to_pfam.keys()]
#Retain phosphorylation interactions
kinases = gname_to_pfam.keys()
BioGRID_A = BioGRID_subset[BioGRID_subset['Official Symbol Interactor A'].isin(kinases)]
BioGRID_B = BioGRID_subset[BioGRID_subset['Official Symbol Interactor B'].isin(kinases)]
BioGRID_subset = pd.concat([BioGRID_A, BioGRID_B])
# region kinase-kinase interactions
# #Remove kinase-kinase interactions
# BioGRID_format = BioGRID_subset[~((BioGRID_subset['Official Symbol Interactor A'].isin(kinases)) & (BioGRID_subset['Official Symbol Interactor B'].isin(kinases)))]
# endregion
BioGRID_format = BioGRID_subset.copy()
BioGRID_format.drop_duplicates(inplace=True)
#There are missing UniProt IDs, but there are no missing gene names.
# print(BioGRID_format.apply(lambda x: (x == '-').sum()))
#Map gene names to UniProt IDs
#Convert all the UniProt IDs as lists
BioGRID_format["protein1"] = BioGRID_format["protein1"].apply(lambda x: x.split(";") if x!='-' else x)
BioGRID_format["protein2"] = BioGRID_format["protein2"].apply(lambda x: x.split(";") if x!='-' else x)
# region gname_to_uniprot
# genes_A = BioGRID_format.loc[BioGRID_format['protein1'] == '-', 'Official Symbol Interactor A'].unique()
# genes_B = BioGRID_format.loc[BioGRID_format['protein2'] == '-', 'Official Symbol Interactor B'].unique()
# gene_names = list(set(list(genes_A) + list(genes_B)))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# # gname_to_uniprot = {}
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=accession'
#     response = requests.get(query_url)
#     uniprot = response.text.split('\n')[1:-1]
#     gname_to_uniprot[gene_name] = uniprot
# endregion
BioGRID_format['protein1'] = BioGRID_format.apply(lambda x: tuple(gname_to_uniprot[x['Official Symbol Interactor A']]) if x['protein1'] == '-' else tuple(x['protein1']), axis=1) #To overcome TypeError: unhashable type: 'list' while dropping duplicates
BioGRID_format['protein2'] = BioGRID_format.apply(lambda x: tuple(gname_to_uniprot[x['Official Symbol Interactor B']]) if x['protein2'] == '-' else tuple(x['protein2']), axis=1)
#Drop the rows that do not have UniProt IDs (pseudogenes or RNA genes)
BioGRID_format = BioGRID_format[~((BioGRID_format["protein1"].apply(len) == 0) | (BioGRID_format["protein2"].apply(len) == 0))]
#BioGRID contains bidirectional interactions. Drop them to get exact kinase-substrate interactions
BioGRID_format = BioGRID_format[~BioGRID_format.apply(lambda x: tuple(sorted(x[['protein1', 'protein2']])), axis=1).duplicated(keep='first')]
# Check for autophosphorylation
# kk = BioGRID_format[BioGRID_format['protein1']==BioGRID_format['protein2']] # 276 interactions
# with open('.\Variables\BioGRID_format.pkl', 'wb') as f:
#     pickle.dump(BioGRID_format, f)

#Overlaps
#BioGRID vs prediction tools
merge1 = pd.concat(BioGRID_format.apply(check_condition, args=(overlaps,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(overlaps.columns) + list(BioGRID_format.columns))
merge2 = pd.concat(BioGRID_format.apply(check_condition_reverse, args=(overlaps,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(overlaps.columns) + list(BioGRID_format.columns))
BioGRID_overlaps = pd.concat([merge1, merge2])
BioGRID_overlaps = BioGRID_overlaps.drop(['protein1', 'protein2','KinaseID_y'], axis=1)
#No duplicates to drop
BioGRID_overlaps.to_csv(r'..\Results\BioGRIDvsPredictiontools.tsv', sep="\t", index=False)

#BioGRID vs NetworKIN and Phospho.ELM database
merge1 = pd.concat(BioGRID_format.apply(check_condition, args=(overlaps_NetworKIN,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(overlaps_NetworKIN.columns) + list(BioGRID_format.columns))
merge2 = pd.concat(BioGRID_format.apply(check_condition_reverse, args=(overlaps_NetworKIN,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(overlaps_NetworKIN.columns) + list(BioGRID_format.columns))
BioGRID_NetworKIN = pd.concat([merge1, merge2])
BioGRID_NetworKIN = BioGRID_NetworKIN.drop(['protein1', 'protein2'], axis=1)
#No duplicates to drop
BioGRID_NetworKIN.to_csv(r'..\Results\BioGRIDvsOverlapsNetworKIN.tsv', sep="\t", index=False)

#BioGRID vs GPS and PhosphoSitePlus database
merge1 = pd.concat(BioGRID_format.apply(check_condition, args=(overlaps_GPS,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(overlaps_GPS.columns) + list(BioGRID_format.columns))
merge2 = pd.concat(BioGRID_format.apply(check_condition_reverse, args=(overlaps_GPS,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(overlaps_GPS.columns) + list(BioGRID_format.columns))
BioGRID_GPS = pd.concat([merge1, merge2])
BioGRID_GPS = BioGRID_GPS.drop(['protein1', 'protein2'], axis=1)
#No duplicates to drop
BioGRID_GPS.to_csv(r'..\Results\BioGRIDvsOverlapsGPS.tsv', sep="\t", index=False)

#Pathway Commons
PC = pd.read_csv(r'..\Data\PathwayCommons12.All.hgnc.txt', sep="\t")
#Retain phosphorylation interactions
# print(PC['INTERACTION_TYPE'].unique())
PC_subset = PC[PC['INTERACTION_TYPE'].isin(['controls-state-change-of',
                                           'controls-transport-of',
                                           'controls-phosphorylation-of',
                                           'controls-expression-of',
                                           'interacts-with'])]
#Map gene names to protein families to filter out those which have "kinase"
# region gname_to_pfam
# genes= list(set(PC_subset['PARTICIPANT_A'].tolist()+PC_subset['PARTICIPANT_B'].tolist()))
# gene_names = [x for x in genes if (x not in gname_to_pfam.keys()) and (x not in non_kinase_genes)]
# base_url = 'https://rest.uniprot.org/uniprotkb'
# for gene_name in gene_names[1000:]:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url)
#     protein_name = response.text.split('\n')[1]
#     if "kinase" in protein_name.lower():
#         gname_to_pfam[gene_name] = protein_name
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
# non_kinase_genes.extend([x for x in gene_names if x not in gname_to_pfam.keys()])
# with open(r'.\Variables\non_kinase_genes.pkl', 'wb') as f:
#     pickle.dump(non_kinase_genes, f)
#Retain only the entries that have kinase
kinases = gname_to_pfam.keys()
PC_A = PC_subset[PC_subset['PARTICIPANT_A'].isin(kinases)]
PC_B = PC_subset[PC_subset['PARTICIPANT_B'].isin(kinases)]
PC_subset = pd.concat([PC_A, PC_B])
# region kinase-kinase interactions
# #Remove kinase-kinase interactions
# PC_subset = PC_subset[~((PC_subset['PARTICIPANT_A'].isin(kinases)) & (PC_subset['PARTICIPANT_B'].isin(kinases)))]
# endregion
PC_subset.drop_duplicates(inplace=True)
#Map gene names to UniProt IDs
# region gname_to_uniprot
# genes= list(set(PC_subset['PARTICIPANT_A'].tolist()+PC_subset['PARTICIPANT_B'].tolist()))
# gene_names = [x for x in genes if x not in gname_to_uniprot.keys()]
# base_url = 'https://rest.uniprot.org/uniprotkb'
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=accession'
#     response = requests.get(query_url)
#     uniprot = response.text.split('\n')[1:-1]
#     gname_to_uniprot[gene_name] = uniprot
# with open('.\Variables\gname_to_uniprot.pkl', 'wb') as f:
#     pickle.dump(gname_to_uniprot, f)
# endregion
PC_format = PC_subset.iloc[:,[0,2,3,1]]
PC_format = PC_format.copy() #To overcome "SettingWithCopyWarning" message
PC_format[['protein1', 'protein2']] = PC_format[['PARTICIPANT_A', 'PARTICIPANT_B']].applymap(lambda x: tuple(gname_to_uniprot[x])) #To overcome TypeError: unhashable type: 'list' while dropping duplicates
PC_format = PC_format.iloc[:,[0,1,4,5,2,3]]
#Drop the rows that do not have UniProt IDs (pseudogenes or RNA genes)
PC_format = PC_format[~((PC_format["protein1"].apply(len) == 0) | (PC_format["protein2"].apply(len) == 0))]
#The interactions belong to only 5 of the 14 binary interaction types
# print(PC_format['INTERACTION_TYPE'].unique())
# with open(r'.\Variables\PC_format.pkl', 'wb') as f:
#     pickle.dump(PC_format, f)

#Overlaps
#PC vs prediction tools
merge1 = pd.concat(PC_format.apply(check_condition, args=(overlaps,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(overlaps.columns) + list(PC_format.columns))
merge2 = pd.concat(PC_format.apply(check_condition_reverse, args=(overlaps,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(overlaps.columns) + list(PC_format.columns))
PC_overlaps = pd.concat([merge1, merge2])
PC_overlaps = PC_overlaps.drop(['protein1', 'protein2','KinaseID_y'], axis=1)
#No duplicates to drop
PC_overlaps.to_csv(r'..\Results\PCvsPredictiontools.tsv', sep="\t", index=False)

#PC vs NetworKIN and Phospho.ELM database
merge1 = pd.concat(PC_format.apply(check_condition, args=(overlaps_NetworKIN,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(overlaps_NetworKIN.columns) + list(PC_format.columns))
merge2 = pd.concat(PC_format.apply(check_condition_reverse, args=(overlaps_NetworKIN,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(overlaps_NetworKIN.columns) + list(PC_format.columns))
PC_NetworKIN = pd.concat([merge1, merge2])
PC_NetworKIN = PC_NetworKIN.drop(['protein1', 'protein2'], axis=1)
#No duplicates to drop
PC_NetworKIN.to_csv(r'..\Results\PCvsOverlapsNetworKIN.tsv', sep="\t", index=False)

#PC vs GPS and PhosphoSitePlus database
merge1 = pd.concat(PC_format.apply(check_condition, args=(overlaps_GPS,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(overlaps_GPS.columns) + list(PC_format.columns))
merge2 = pd.concat(PC_format.apply(check_condition_reverse, args=(overlaps_GPS,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(overlaps_GPS.columns) + list(PC_format.columns))
PC_GPS = pd.concat([merge1, merge2])
PC_GPS = PC_GPS.drop(['protein1', 'protein2'], axis=1)
#No duplicates
PC_GPS.to_csv(r'..\Results\PCvsOverlapsGPS.tsv', sep="\t", index=False)

##Number of vertices and edges
#Remove bidirectional interactions to get exact kinase-substrate interactions
PC_total = PC_format[~PC_format.apply(lambda x: tuple(sorted(x[['protein1', 'protein2']])), axis=1).duplicated(keep='first')]
# Check for autophosphorylation
# kk = PC_total[PC_total['protein1']==PC_total['protein2']] # No interactions

net_prop = pd.concat([net_prop,
                     pd.DataFrame({'Sources': ['STRING', 'BioGRID', 'Pathway Commons'],
                                   'Vertices': [pd.concat([STRING_format['protein1'],STRING_format['protein2']]).nunique(),
                                                pd.concat([BioGRID_format['protein1'], BioGRID_format['protein2']]).nunique(),
                                                pd.concat([PC_total['protein1'], PC_total['protein2']]).nunique()],
                                   'Edges': [len(STRING_format.index),len(BioGRID_format.index),len(PC_total.index)]}).set_index("Sources")])
net_prop = net_prop[~net_prop.index.duplicated(keep='last')]
# net_prop.to_csv(r'..\Results\network_properties.tsv', sep="\t")

##Venn diagrams for the visualisations of the overlaps
#Labels are the number of edges (predictions)
matplotlib.use('Qt5Agg') #To prevent matplotlib from not responding
#STRING vs prediction tools
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(STRING_format.index)*0.01, len(overlaps.index), len(STRING_overlaps.index)),
             set_colors = ('green', 'orange'),
             set_labels = ('STRING', 'Prediction tools'))
ax.set_title('STRING vs prediction tools \n (GPS and NetworKIN)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['STRING','Edges']-len(STRING_overlaps.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps.index)-len(STRING_overlaps.index)}")
venn.get_label_by_id('11').set_text(f"{len(STRING_overlaps.index)}")
venn.get_label_by_id('11').set_position((0.43,0))
venn.get_label_by_id('01').set_position((0.59,0))
venn.set_labels[1].set_x(0.53)
venn.set_labels[1].set_y(-0.18)
plt.legend([len(STRING_format.index), len(overlaps.index)],title="Total interactions",bbox_to_anchor=(0.85, 1.01))
plt.savefig(r'..\Results\Figures\STRINGvsPredictiontools.png')

#STRING vs NetworKIN and Phospho.ELM database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(STRING_format.index)*0.01, len(overlaps_NetworKIN.index), len(STRING_NetworKIN.index)),
             set_colors = ('green', 'skyblue'),
             set_labels = ('STRING', 'NetworKIN ∩ Phospho.ELM'))
ax.set_title('STRING vs NetworKIN and Phospho.ELM database')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['STRING','Edges']-len(STRING_NetworKIN.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps_NetworKIN.index)-len(STRING_NetworKIN.index)}")
venn.get_label_by_id('11').set_text(f"{len(STRING_NetworKIN.index)}")
venn.get_label_by_id('11').set_position((0.49,0))
venn.get_label_by_id('01').set_position((0.59,0))
venn.set_labels[1].set_x(0.52)
plt.legend([len(STRING_format.index), len(overlaps_NetworKIN.index)],title="Total interactions",bbox_to_anchor=(0.9, 0.95))
plt.savefig('..\Results\Figures\STRINGvsOverlapsNetworKIN.png')

#STRING vs GPS and PhosphoSitePlus database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(STRING_format.index)*0.1, len(overlaps_GPS.index), len(STRING_GPS.index)),
             set_colors = ('green', 'purple'),
             set_labels = ('STRING', 'GPS ∩ PhosphoSitePlus'))
ax.set_title('STRING vs GPS and PhosphoSitePlus database')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['STRING','Edges']-len(STRING_GPS.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps_GPS.index)-len(STRING_GPS.index)}")
venn.get_label_by_id('11').set_text(f"{len(STRING_GPS.index)}")
plt.legend([len(STRING_format.index), len(overlaps_GPS.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.98))
venn.set_labels[1].set_x(0.39)
plt.savefig(r'..\Results\Figures\STRINGvsOverlapsGPS.png')

#BioGRID vs prediction tools
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(BioGRID_format.index)*0.1, len(overlaps.index), len(BioGRID_overlaps.index)),
             set_colors = ('red', 'orange'),
             set_labels = ('BioGRID', 'Prediction tools'))
ax.set_title('BioGRID vs prediction tools \n (GPS and NetworKIN)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['BioGRID','Edges']-len(BioGRID_overlaps.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps.index)-len(BioGRID_overlaps.index)}")
venn.get_label_by_id('11').set_text(f"{len(BioGRID_overlaps.index)}")
venn.set_labels[1].set_x(0.56)
venn.set_labels[1].set_y(-0.15)
plt.legend([len(BioGRID_format.index), len(overlaps.index)],title="Total interactions",bbox_to_anchor=(0.85, 1))
plt.savefig(r'..\Results\Figures\BioGRIDvsPredictiontools.png')

#BioGRID vs NetworKIN and Phospho.ELM database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(BioGRID_format.index)*0.1, len(overlaps_NetworKIN.index), len(BioGRID_NetworKIN.index)),
             set_colors = ('red', 'skyblue'),
             set_labels = ('BioGRID', 'NetworKIN ∩ Phospho.ELM'))
ax.set_title('BioGRID vs NetworKIN and Phospho.ELM database')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['BioGRID','Edges']-len(BioGRID_NetworKIN.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps_NetworKIN.index)-len(BioGRID_NetworKIN.index)}")
venn.get_label_by_id('11').set_text(f"{len(BioGRID_NetworKIN.index)}")
venn.set_labels[1].set_x(0.52)
plt.legend([len(BioGRID_format.index), len(overlaps_NetworKIN.index)],title="Total interactions",bbox_to_anchor=(0.9, 0.95))
plt.savefig(r'..\Results\Figures\BioGRIDvsOverlapsNetworKIN.png')

#BioGRID vs GPS and PhosphoSitePlus database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(BioGRID_format.index), len(overlaps_GPS.index), len(BioGRID_GPS.index)),
             set_colors = ('red', 'purple'),
             set_labels = ('BioGRID', 'GPS ∩ PhosphoSitePlus'))
ax.set_title('BioGRID vs GPS and PhosphoSitePlus database')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['BioGRID','Edges']-len(BioGRID_GPS.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps_GPS.index)-len(BioGRID_GPS.index)}")
venn.get_label_by_id('11').set_text(f"{len(BioGRID_GPS.index)}")
venn.set_labels[1].set_x(0.48)
plt.legend([len(BioGRID_format.index), len(overlaps_GPS.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.98))
plt.savefig(r'..\Results\Figures\BioGRIDvsOverlapsGPS.png')

#Pathway Commons vs prediction tools
PC_overlaps = PC_overlaps.drop(['INTERACTION_DATA_SOURCE','INTERACTION_TYPE','PARTICIPANT_A','PARTICIPANT_B'], axis=1)
PC_overlaps.drop_duplicates(inplace=True)
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PC_total.index)*0.1, len(overlaps.index), len(PC_overlaps.index)),
             set_colors = ('grey', 'orange'),
             set_labels = ('Pathway Commons', 'Prediction tools'))
ax.set_title('Pathway Commons vs prediction tools \n (GPS and NetworKIN)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Pathway Commons','Edges']-len(PC_overlaps.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps.index)-len(PC_overlaps.index)}")
venn.get_label_by_id('11').set_text(f"{len(PC_overlaps.index)}")
venn.set_labels[1].set_x(0.56)
venn.set_labels[1].set_y(-0.15)
plt.legend([len(PC_total.index), len(overlaps.index)],title="Total interactions",bbox_to_anchor=(0.85, 1))
plt.savefig(r'..\Results\Figures\PCvsPredictiontools.png')

#Pathway Commons vs NetworKIN and Phospho.ELM database
PC_NetworKIN = PC_NetworKIN.drop(['INTERACTION_DATA_SOURCE','INTERACTION_TYPE','PARTICIPANT_A','PARTICIPANT_B'], axis=1)
PC_NetworKIN.drop_duplicates(inplace=True)
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PC_total.index)*0.1, len(overlaps_NetworKIN.index), len(PC_NetworKIN.index)),
             set_colors = ('grey', 'skyblue'),
             set_labels = ('Pathway Commons', 'NetworKIN ∩ Phospho.ELM'))
ax.set_title('Pathway Commons vs NetworKIN and Phospho.ELM database')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Pathway Commons','Edges']-len(PC_NetworKIN.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps_NetworKIN.index)-len(PC_NetworKIN.index)}")
venn.get_label_by_id('11').set_text(f"{len(PC_NetworKIN.index)}")
venn.set_labels[1].set_x(0.52)
plt.legend([len(PC_total.index), len(overlaps_NetworKIN.index)],title="Total interactions",bbox_to_anchor=(0.9, 0.95))
plt.savefig(r'..\Results\Figures\PCvsOverlapsNetworKIN.png')

#Pathway Commons vs GPS and PhosphoSitePlus database
PC_GPS = PC_GPS.drop(['INTERACTION_DATA_SOURCE','INTERACTION_TYPE','PARTICIPANT_A','PARTICIPANT_B'], axis=1)
PC_GPS.drop_duplicates(inplace=True)
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PC_total.index), len(overlaps_GPS.index), len(PC_GPS.index)),
             set_colors = ('grey', 'purple'),
             set_labels = ('Pathway Commons', 'GPS ∩ PhosphoSitePlus'))
ax.set_title('Pathway Commons vs GPS and PhosphoSitePlus database')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Pathway Commons','Edges']-len(PC_GPS.index)}")
venn.get_label_by_id('01').set_text(f"{len(overlaps_GPS.index)-len(PC_GPS.index)}")
venn.get_label_by_id('11').set_text(f"{len(PC_GPS.index)}")
venn.set_labels[1].set_x(0.48)
plt.legend([len(PC_total.index), len(overlaps_GPS.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.98))
plt.savefig(r'..\Results\Figures\PCvsOverlapsGPS.png')