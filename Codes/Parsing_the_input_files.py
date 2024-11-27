import pandas as pd
import numpy as np
from unipressed import IdMappingClient, UniprotkbClient
import time
import requests
import copy
import pickle
from modules import find_overlaps, find_overlaps2
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import matplotlib

#Load saved variables
with open('.\Variables\gname_to_uniprot.pkl', 'rb') as f:
    gname_to_uniprot = pickle.load(f)
with open('.\Variables\ensp_to_uniprot.pkl', 'rb') as f:
    ensp_to_uniprot = pickle.load(f)
with open('.\Variables\hproteins.pkl', 'rb') as f:
    hproteins = pickle.load(f)
with open(r'.\Variables\uniprot_to_pfam.pkl', 'rb') as f:
    uniprot_to_pfam = pickle.load(f)
with open(r'.\Variables\gname_to_pfam.pkl', 'rb') as f:
    gname_to_pfam = pickle.load(f)

##Read the xlsx files of the predictions
GPS = pd.read_excel(r"..\Data\GPS_predictions.xlsx", sheet_name=1) #sheet_name=1 reads the second sheet in the excel file
# print(GPS.head(5))
NetworKIN = pd.read_excel(r"..\Data\NetworKIN_predictions.xlsx")
# print(NetworKIN.head(5))

#For GPS, extract information for only homo sapiens
GPS_subset = GPS[GPS['Species']=='Homo sapiens']
# Select only the substrates that belong to humans
# region human_proteins
# ids = GPS_subset['UniProt ID (Substrate)'].unique().tolist()
# base_url = 'https://rest.uniprot.org/uniprotkb'
# # hproteins=[]
# for id in ids[3000:]:
#     query_url = f'{base_url}/search?format=tsv&query=accession:{id}&fields=organism_id'
#     response = requests.get(query_url, timeout=10)
#     org = response.text.split('\n')[1]
#     if org=="9606":
#         hproteins.append(id)
# with open(r'.\Variables\hproteins.pkl', 'wb') as f:
#     pickle.dump(hproteins, f)
# endregion
GPS_subset = GPS_subset[GPS_subset['UniProt ID (Substrate)'].isin(hproteins)]
# print(GPS_subset.isna().sum())
#Only UniProt IDs of kinases have missing values.
#print(GPS_subset.head(5))

##Parse the files to columns of Uniprot ID, Position, Residue, Kinase and Peptide sequence
#GPS
GPS_format = GPS_subset.iloc[:, [2,3,4,0,5,1]]
GPS_format.columns=['SubstrateID','Position','Code','KinaseID','Peptide']+GPS_format.columns[5:].tolist()
GPS_format = GPS_format.copy() #To overcome "SettingWithCopyWarning" message
GPS_format['KinaseID'] = GPS_format['KinaseID'].apply(lambda x: x.split(";") if str(x)!='nan' else x)
#Map gene names to UniProt IDs
# region gname_to_uniprot
# genes= GPS_format['Gene name'][GPS_format['UniProt ID (Kinase)'].isna()].unique().tolist()
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
GPS_format['KinaseID'] = GPS_format.apply(lambda x: tuple(gname_to_uniprot[x['Gene name']]) if str(x['KinaseID'])=='nan' else tuple(x['KinaseID']), axis=1) #To overcome TypeError: unhashable type: 'list' while dropping duplicates
GPS_format =GPS_format[GPS_format['KinaseID']!=()]
GPS_format.loc[:,'Peptide']=GPS_format['Peptide'].apply(lambda x: x[:30]+x[30].lower()+x[31:]) #Represent the phosphorylation site in lowercase
# GPS_format['Peptide'].str.len().mean()
GPS_format.loc[:,'Peptide']=GPS_format['Peptide'].str.slice(start=26, stop=35) #Convert the 61-mer to 9-mer
#I validated this representation by comparing the peptide sequence of an overlap between GPS and NetworKIN (highlighted in green in Excel)
# For positions less than 5, retain only the right end of the phosphorylation site (5-mer)
GPS_format.loc[GPS_format['Position']<5, 'Peptide']=GPS_format['Peptide'].str.slice(start=4, stop=9)
#If there is '*' in the 9-mer, slice it out
GPS_format['Peptide'] = GPS_format['Peptide'].str.replace('*', '', regex=True)
# region kinase-kinase interactions
# #Extract only kinase-substrate interactions
# #Map uniprot ids of substrates to protein familes to filter out those which have "kinase"
# # region uniprot_to_pfam
# # uniprot= GPS_format['SubstrateID'].unique().tolist()
# # base_url = 'https://rest.uniprot.org/uniprotkb'
# # # uniprot_to_pfam={}
# # for id in uniprot[2000:]:
# #     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
# #     response = requests.get(query_url, timeout=10)
# #     protein_fam = response.text.split('\n')[1]
# #     if "kinase" in protein_fam.lower():
# #         uniprot_to_pfam[id] = protein_fam
# # with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
# #     pickle.dump(uniprot_to_pfam, f)
# # endregion
# kinases = uniprot_to_pfam.keys()
# #Remove kinase-kinase interactions
# GPS_format = GPS_format[~GPS_format['SubstrateID'].isin(kinases)]
# endregion
#Check for proteins that are not kinases in the kinase column
# region gname_to_pfam
# gene_names = GPS_format['Gene name'].unique().tolist()
# base_url = 'https://rest.uniprot.org/uniprotkb'
# # gname_to_pfam={}
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     if "kinase" in protein_fam.lower():
#         gname_to_pfam[gene_name] = protein_fam
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
# gene_names = [x for x in gene_names if x not in gname_to_pfam.keys()]
#8 proteins that are not kinases were found
#They were verified individually and added to the dictionary
# region gname_to_pfam
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     gname_to_pfam[gene_name] = protein_fam
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
# Check for autophosphorylation
# GPS_format['flag'] = GPS_format.apply(lambda row: find_overlaps2(row['SubstrateID'], row['KinaseID']), axis=1)
# kk = GPS_format.dropna(subset="flag") # 852 interactions
GPS_format.drop_duplicates(inplace=True)
# GPS_format.to_csv(r'.\Variables\GPS_format.tsv', sep="\t", index=False)

#NetworKIN
#Convert Ensembl protein IDs of kinase to Uniprot IDs
# region ensp_to_uniprot
# #Add all the entries of BioMart to 'ensp_to_uniprot' file
# biomart = pd.read_csv('..\Data\mart_export.txt', sep='\t')
# biomart.dropna(inplace=True)
# biomart.drop_duplicates(inplace=True)
# #There are 14 ENSPs that have more than one corresponding UniProt IDs.
# duplicate = list(set(biomart['Protein stable ID'][biomart.duplicated(['Protein stable ID'])]))
# #Drop the above 14 ENSP entries (46 in total) and save the rest as a dictionary
# biomart_unique = biomart[~biomart['Protein stable ID'].isin(duplicate)]
# ensp_to_uniprot = dict(zip(biomart_unique['Protein stable ID'], biomart_unique['UniProtKB/Swiss-Prot ID']))
# #Convert the 14 ENSPs using unipressed
# request = IdMappingClient.submit(
#     source="Ensembl_Protein", dest="UniProtKB", ids=duplicate
# )
# time.sleep(5.0)
# converted = list(request.each_result())
# ensp = list(map(lambda x: x['from'], converted))
# uniprot = list(map(lambda x: x['to'], converted))
# for i in range(len(ensp)):
#     ensp_to_uniprot[ensp[i]]= uniprot[i]
#
# #Unipressed
# request = IdMappingClient.submit(
#     source="Ensembl_Protein", dest="UniProtKB", ids=NetworKIN["ensemblID_kinase"].to_list()
# )
# time.sleep(5.0)
# # print(list(request.each_result()))
# #NetworKIN["ensemblID_kinase"].isnull().sum()
# #No empty cells.
# # print(len(list(request.each_result())))
# converted = list(request.each_result())
# ensp = list(map(lambda x: x['from'], converted))
# uniprot = list(map(lambda x: x['to'], converted))
# for i in range(len(ensp)):
#     ensp_to_uniprot[ensp[i]]= uniprot[i]
# endregion
NetworKIN_format= copy.deepcopy(NetworKIN)
NetworKIN_format['ensemblID_kinase'] = NetworKIN_format['ensemblID_kinase'].replace(ensp_to_uniprot)
#Predictions whose Ensembl Protein IDs of kinase has not been converted to UniProt IDs, are assigned with 'nan'.
NetworKIN_format['ensemblID_kinase'] = NetworKIN_format['ensemblID_kinase'].apply(lambda x: np.nan if x.startswith('EN') else x)
# NetworKIN_format.isna().sum()
#The conversion to UniProt IDs was done for only 4038 out of 7143 entries.
#Map gene names to UniProt IDs
NetworKIN_format['ensemblID_kinase'] = NetworKIN_format['ensemblID_kinase'].apply(lambda x: x.split(";") if str(x)!='nan' else x)
# region gname_to_uniprot
# genes= NetworKIN_format['genesymbol_kinase'][NetworKIN_format['ensemblID_kinase'].isna()].unique().tolist()
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
NetworKIN_format.dropna(subset=["ensemblID_kinase", "genesymbol_kinase"], how="all", inplace=True) #'all' indicates that only rows where both columns are nan should be dropped.
NetworKIN_format['ensemblID_kinase'] = NetworKIN_format.apply(lambda x: tuple(gname_to_uniprot[x['genesymbol_kinase']]) if str(x['ensemblID_kinase'])=='nan' else tuple(x['ensemblID_kinase']), axis=1) #To overcome TypeError: unhashable type: 'list' while dropping duplicates
# NetworKIN_format = NetworKIN_format[NetworKIN_format['ensemblID_kinase']!=()]
NetworKIN_format['Code'] = NetworKIN_format["residue"].apply(lambda x: str(x)[0])
NetworKIN_format['Position'] = NetworKIN_format["residue"].apply(lambda x: int(str(x)[1:]))
NetworKIN_format = NetworKIN_format.iloc[:, [0,15,14,6,11,8,3,4]] #addition of motif and context scores
NetworKIN_format.columns=['SubstrateID','Position','Code','KinaseID','Peptide']+NetworKIN_format.columns[5:].tolist()
# region kinase-kinase interactions
# #Extract only kinase-substrate interactions
# #Map uniprot ids of substrates to protein familes to filter out those which have "kinase"
# # region uniprot_to_pfam
# # uniprot= NetworKIN_format['SubstrateID'].unique().tolist()
# # base_url = 'https://rest.uniprot.org/uniprotkb'
# # uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# # for id in uniprot[1000:]:
# #     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
# #     response = requests.get(query_url, timeout=10)
# #     protein_fam = response.text.split('\n')[1]
# #     if "kinase" in protein_fam.lower():
# #         uniprot_to_pfam[id] = protein_fam
# # with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
# #     pickle.dump(uniprot_to_pfam, f)
# # endregion
# kinases = uniprot_to_pfam.keys()
# #Remove kinase-kinase interactions
# NetworKIN_format = NetworKIN_format[~NetworKIN_format['SubstrateID'].isin(kinases)]
# endregion
#Check for proteins that are not kinases in the kinase column
# region gname_to_pfam
# gene_names = NetworKIN_format['genesymbol_kinase'].unique().tolist()
# base_url = 'https://rest.uniprot.org/uniprotkb'
# gene_names = [x for x in gene_names if x not in gname_to_pfam.keys()]
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     if "kinase" in protein_fam.lower():
#         gname_to_pfam[gene_name] = protein_fam
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
# gene_names = [x for x in gene_names if x not in gname_to_pfam.keys()]
# print(len(NetworKIN_format[NetworKIN_format['genesymbol_kinase'].isna()]))
#89 interactions have missing gene names, but the uniprot IDs of kinases are available.
uniprot = NetworKIN_format.loc[NetworKIN_format['genesymbol_kinase'].isna(), 'KinaseID'].unique().tolist()
uniprot = [str(x[0]) for x in uniprot] # Convert tuple of strings to strings
uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# All are protein kinases in the kinase column
# Check for autophosphorylation
# NetworKIN_format['flag'] = NetworKIN_format.apply(lambda row: find_overlaps2(row['SubstrateID'], row['KinaseID']), axis=1)
# kk = NetworKIN_format.dropna(subset="flag") # No interactions
# NetworKIN_format.to_csv(r'.\Variables\NetworKIN_format.tsv', sep="\t", index=False)

##Overlaps of predictions
#Ignoring positions
s1 = pd.merge(NetworKIN_format, GPS_format, how='inner', on=['SubstrateID','Code','Peptide'])
#Positions represent the correct phosphorylation sites. But the position numbers are slightly different in some of the predictions.
s1['KinaseID'] = s1.apply(lambda row: find_overlaps(row['KinaseID_x'], row['KinaseID_y']), axis=1)
s1 = s1[s1['KinaseID'].apply(lambda x: len(x) > 0)]
s1.sort_values(by=['motif_score', 'context_score'], inplace=True, ascending=False)
#The motif scores are above the threshold of 0.5. Filter out the predictions with context scores <0.9.
#Three of the hits has a large difference in the position numbers, so limit the difference to be 2
s1['filter'] = s1.apply(lambda x: (abs(x['Position_x']-x['Position_y'])<=2) and (x['context_score']>=0.9), axis=1)
# s1.sort_values(by=['Peptide'], key=lambda x: x.str.len(), inplace=True)
overlaps = s1[s1['filter'] == True].iloc[: , :-1]
order = [0,11,2,1,4,8,6,7,3,5,9,10]
overlaps = overlaps[[overlaps.columns[i] for i in order]]
overlaps.to_csv(r'..\Results\overlaps_phosphorylation.tsv', sep="\t", index=False)

# region Code testing
# s1 = pd.merge(NetworKIN_format, NetworKIN_format, how='inner', on=['SubstrateID','Code','Peptide']) # (or)
# s1 = pd.merge(GPS_format, GPS_format, how='inner', on=['SubstrateID','Code','Peptide'])
# s1['KinaseID'] = s1.apply(lambda row: find_overlaps(row['KinaseID_x'], row['KinaseID_y']), axis=1)
# s1 = s1[s1['KinaseID'].apply(lambda x: len(x) > 0)]
# s1['filter'] = s1.apply(lambda x: (abs(x['Position_x']-x['Position_y'])<=2), axis=1)
# overlaps = s1[s1['filter'] == True].iloc[: , :-1]
# endregion

#Perfect match
s2 = copy.deepcopy(overlaps)
s2['filter'] = s2.apply(lambda x: x['Position_x']==x['Position_y'], axis=1)
s2 = s2[s2['filter'] == True].iloc[: , :-1]
# region Code testing (Cont.)
# order = [0,1,2,9,4,5] #Code testing (GPS_format)
# s2 = s2[[s2.columns[i] for i in order]]
# s2 = s2.rename(columns={s2.columns[1]: "Position", s2.columns[5]: "Gene name"})
# s2.drop_duplicates(inplace=True)
# endregion
order = [0,3,2,1,4,6,7,8,9,10,11]
s2 = s2[[s2.columns[i] for i in order]]
s2 = s2.rename(columns={s2.columns[1]: "Position"})
s2.to_csv(r'..\Results\perfect_match_phosphorylation.tsv', sep="\t", index=False)

##Check for overlaps between the predictions and the databases used in the papers
#NetworKIN and Phospho.ELM database
PhosphoELM = pd.read_csv("..\Data\phosphoELM_vertebrate_2015-04.dump", sep="\t")
PhosphoELM_subset = PhosphoELM[PhosphoELM['species']=='Homo sapiens'] #Extract information for only homo sapiens
# region human_proteins
# ids = PhosphoELM_subset['acc'].unique().tolist()
# base_url = 'https://rest.uniprot.org/uniprotkb'
# ids = [x for x in ids if x not in hproteins]
# for id in ids[3000:]:
#     query_url = f'{base_url}/search?format=tsv&query=accession:{id}&fields=organism_id'
#     response = requests.get(query_url, timeout=10)
#     org = response.text.split('\n')[1]
#     if org=="9606":
#         hproteins.append(id)
# with open(r'.\Variables\hproteins.pkl', 'wb') as f:
#     pickle.dump(hproteins, f)
# endregion
PhosphoELM_subset = PhosphoELM_subset[PhosphoELM_subset['acc'].isin(hproteins)]

#Convert Ensembl protein IDs of substrate to UniProt IDs
# region ensp_to_uniprot
# ensp = PhosphoELM_subset['acc'].loc[PhosphoELM_subset['acc'].str.startswith('EN')].tolist()
# request = IdMappingClient.submit(
#     source="Ensembl_Protein", dest="UniProtKB", ids=ensp
# )
# time.sleep(5.0)
# converted = list(request.each_result())
# # ensp = set(ensp)
# #There are 331 unique Ensembl protein IDs, out of which only 245 were converted to UniProt IDs.
# ensp = list(map(lambda x: x['from'], converted))
# uniprot = list(map(lambda x: x['to'], converted))
# for i in range(len(ensp)):
#     ensp_to_uniprot[ensp[i]]= uniprot[i]
# with open('.\Variables\ensp_to_uniprot.pkl', 'wb') as f:
#     pickle.dump(ensp_to_uniprot, f)
# endregion
PhosphoELM_subset= PhosphoELM_subset.copy() #To overcome "SettingWithCopyWarning" message
PhosphoELM_subset['acc'] = PhosphoELM_subset['acc'].replace(ensp_to_uniprot)
#UniProt IDs of kinases obtained via web scrapping
tables_on_page = pd.read_html("..\Data\Phospho.ELM Result.html")
table = tables_on_page[1]
#Remove the predictions whose Ensembl Protein IDs of substrate has not been converted to UniProt IDs.
PhosphoELM_format = PhosphoELM_subset.drop(PhosphoELM_subset[PhosphoELM_subset['acc'].str.startswith('EN')].index)
#No interactions are dropped
PhosphoELM_format['sequence']=PhosphoELM_format.apply(lambda x: x['sequence'][:x['position']-1]+x['sequence'][x['position']-1:x['position']].lower()+x['sequence'][x['position']:], axis=1) #Represent the phosphorylation site in lowercase
#'axis=1' is used to apply the lambda function row-wise, without which, the function would be applied column-wise.
# For positions less than 5, retain only the right end of the phosphorylation site (5-mer), and convert the rest to 9-mer based on the 'position' value
PhosphoELM_format['sequence'] = PhosphoELM_format.apply(lambda x: str(x['sequence'])[(x.position-5):(x.position+4)] if x.position >= 5 else str(x['sequence'])[(x.position-1):(x.position+4)], axis=1)
#Add the "Uniprot" and "Description" columns from "table"
PhosphoELM_format = pd.merge(PhosphoELM_format, table[['Name', 'Uniprot','Description']], left_on='kinases', right_on='Name', how='left')
PhosphoELM_format.drop('Name', axis=1, inplace=True)
PhosphoELM_format = PhosphoELM_format.iloc[:, [0,2,3,9,1,9,5,10]]
PhosphoELM_format.columns=['SubstrateID','Position','Code','KinaseID','Peptide','KinaseID_y','Kinases']+PhosphoELM_format.columns[7:].tolist()
# region kinase-kinase interactions
# #Extract only kinase-substrate interactions
# #Map uniprot ids of substrates to protein familes to filter out those which have "kinase"
# # region uniprot_to_pfam
# # uniprot= PhosphoELM_format['SubstrateID'].unique().tolist()
# # base_url = 'https://rest.uniprot.org/uniprotkb'
# # uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# # for id in uniprot[4000:]:
# #     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
# #     response = requests.get(query_url, timeout=10)
# #     protein_fam = response.text.split('\n')[1]
# #     if "kinase" in protein_fam.lower():
# #         uniprot_to_pfam[id] = protein_fam
# # with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
# #     pickle.dump(uniprot_to_pfam, f)
# # endregion
# kinases = uniprot_to_pfam.keys()
# #Remove kinase-kinase interactions
# PhosphoELM_format = PhosphoELM_format[~PhosphoELM_format['SubstrateID'].isin(kinases)]
# endregion
#Check for proteins that are not kinases in the kinase column
# region uniprot_to_pfam
# uniprot= PhosphoELM_format['KinaseID'].unique().tolist()
# base_url = 'https://rest.uniprot.org/uniprotkb'
# uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# for id in uniprot:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     if "kinase" in protein_fam.lower():
#         uniprot_to_pfam[id] = protein_fam
# with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_pfam, f)
# endregion
#3 proteins that are not kinases were found
#They were verified individually. 'P14619' and 'Q7Z370' were added to the dictionary
# query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:"Q13976"&fields=protein_families'
# response = requests.get(query_url, timeout=10)
# protein_fam = response.text.split('\n')[1]
# if "kinase" in protein_fam.lower():
#     uniprot_to_pfam['P14619'] = protein_fam
# query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:"Q7Z370"&fields=protein_families'
# response = requests.get(query_url, timeout=10)
# protein_fam = response.text.split('\n')[1]
# if "kinase" in protein_fam.lower():
#     uniprot_to_pfam['Q7Z370'] = protein_fam
# with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_pfam, f)
kinases = uniprot_to_pfam.keys()
PhosphoELM_format = PhosphoELM_format[PhosphoELM_format['KinaseID'].isin(kinases) | pd.isna(PhosphoELM_format['KinaseID'])]
PhosphoELM_format.drop_duplicates(inplace=True)
# PhosphoELM_format.to_csv(r'.\Variables\PhosphoELM_format.tsv', sep="\t", index=False)

#Overlaps
PhosphoELM_format["KinaseID"] = PhosphoELM_format["KinaseID"].apply(lambda x: tuple(x.split(",")) if str(x)!="nan" else x)
NetworKIN_predictor = NetworKIN_format.join(NetworKIN.loc[:,['predictor']]) #Add the predictor column
order = [*range(5)]+[*range(6,9)]+[5]
NetworKIN_predictor = NetworKIN_predictor[[NetworKIN_predictor.columns[i] for i in order]]
#1. Overlaps for the entries having a UniProt ID for kinases
selected = PhosphoELM_format[~PhosphoELM_format['KinaseID'].isna()]
unselected = PhosphoELM_format[PhosphoELM_format['KinaseID'].isna()]
overlaps_NetworKIN = pd.merge(NetworKIN_predictor, selected, how='inner', on=['SubstrateID','KinaseID','Code','Peptide'])
#2. Overlaps for the entries without UniProt IDs for kinases
unselected= unselected.copy() #To overcome "SettingWithCopyWarning" message
unselected.drop('KinaseID_y',axis=1, inplace=True)
overlaps_unselected = pd.merge(NetworKIN_predictor, unselected, how='inner', on=['SubstrateID','Code','Peptide'])
unselected = overlaps_unselected[~overlaps_unselected['Kinases'].isna()]
selected = overlaps_unselected[overlaps_unselected['Kinases'].isna()]
overlaps_NetworKIN.columns=unselected.columns #Rename "KinaseID" as "KinaseID_x"
overlaps_NetworKIN = pd.concat([overlaps_NetworKIN, selected])
#3. Manually inspect the 'genesymbol_kinase','predictor' and 'Kinases' columns of the unselected entries.
#'genesymbol_kinase'
all_kin_dict = {}
for predictor in unselected["genesymbol_kinase"].unique():
    kinases = unselected.loc[unselected["genesymbol_kinase"] == predictor, "Kinases"].unique()
    all_kin_dict[predictor] = list(kinases)
# df = pd.DataFrame({'key' : all_kin_dict.keys() , 'value' : all_kin_dict.values() })
# df.to_csv(r'.\Variables\all_kin_dict.tsv', sep="\t", index=False)
kin_dict = {"PRKCA":["PKC_group"],
            "AKT1":["PKB_group"],
            "GSK3B":["GSK-3_group"],
            "CSNK2A2":["CK2_group"],
            "MAPK14":["MAPK_group"],
            "SGK": ["SGK_group"],
            "CDK2": ["CDK_group"],
            "MAPK9": ["MAPK_group"],
            "CSNK1D": ["CK1_group"],
            "MAPK1":["MAPK_group"],
            "PRKCG": ["PKC_group"],
            "RPS6KB1": ["p70S6K_group", "RSK_group"],
            "PRKCB1": ["PKC_group"],
            "CDC2": ["CDK_group"]}
not_in_kin_dict={}
for key, value in all_kin_dict.items():
    if key not in kin_dict:
        not_in_kin_dict[key]=value
    else:
        if kin_dict.get(key)==value:
            continue
        not_in_kin_dict[key] = [x for x in value if x not in kin_dict.get(key)]
for key in kin_dict:
    kinase_list = kin_dict.get(key)
    selected = unselected[(unselected['genesymbol_kinase'] == key) & (unselected['Kinases'].isin(kinase_list))]
    overlaps_NetworKIN = pd.concat([overlaps_NetworKIN, selected])
mask = ~unselected[unselected.columns].apply(tuple, axis=1).isin(overlaps_NetworKIN[overlaps_NetworKIN.columns].apply(tuple, axis=1))
#tuple() converts each row of both dataframes to a tuple of values, which we can then compare for equality using isin().
unselected = unselected[mask]

#'predictor'
all_kin_dict = {}
for predictor in unselected["predictor"].unique():
    kinases = unselected.loc[unselected["predictor"] == predictor, "Kinases"].unique()
    all_kin_dict[predictor] = list(kinases)
temp = copy.deepcopy(all_kin_dict)
for key, value in all_kin_dict.items():
    if key in not_in_kin_dict:
        if not_in_kin_dict.get(key)==value:
            del temp[key]
        else:
            temp[key] = [x for x in value if x not in not_in_kin_dict.get(key)]
all_kin_dict = temp
# df = pd.DataFrame({'key' : all_kin_dict.keys() , 'value' : all_kin_dict.values() })
# df.to_csv(r'.\Variables\all_kin_dict.tsv', sep="\t", index=False)
new_info = {"PKB":["PKB_group"],
            "PKA":["PKA_group"],
            "PKG": ["PKG/cGK_group"],
            "CaM-II": ["CaM-KII_group"],
            "RSK": ["RSK_group"]}
for key, value in new_info.items():
    if key in kin_dict:
        kin_dict[key].extend(value)
    else:
        kin_dict[key] = value
for key, value in all_kin_dict.items():
    if key in kin_dict:
        if value == kin_dict.get(key):
            continue
        else:
            if key in not_in_kin_dict:
                not_in_kin_dict[key].extend([x for x in value if x not in kin_dict.get(key)])
            else:
                not_in_kin_dict[key]=[x for x in value if x not in kin_dict.get(key)]
    else:
        if key in not_in_kin_dict:
            not_in_kin_dict[key].extend(value)
        else:
            not_in_kin_dict[key] = value
for key in kin_dict:
    kinase_list = kin_dict.get(key)
    selected = unselected[(unselected['predictor'] == key) & (unselected['Kinases'].isin(kinase_list))]
    overlaps_NetworKIN = pd.concat([overlaps_NetworKIN, selected])
mask = ~unselected[unselected.columns].apply(tuple, axis=1).isin(overlaps_NetworKIN[overlaps_NetworKIN.columns].apply(tuple, axis=1))
unselected = unselected[mask]
# with open('.\Variables\kin_dict.pkl', 'wb') as f:
#     pickle.dump(kin_dict, f)
# with open(r'.\Variables\not_in_kin_dict.pkl', 'wb') as f: #r to overcome "OSError: [Errno 22] Invalid argument: '.\\Variables\not_in_kin_dict.pkl'"
#     pickle.dump(not_in_kin_dict, f)

# region Synonyms
# #Get the synonyms of kinases from UniProtKB (Gene Names column)
# base_url = 'https://rest.uniprot.org/uniprotkb'
# # gene_names = unselected["genesymbol_kinase"].unique()
# gene_names = unselected["predictor"].unique()
# gene_synonyms = {}
# temp=pd.DataFrame()
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=gene_names'
#     response = requests.get(query_url)
#     genes = response.text.split('\n')[1].split(" ")
#     if genes==['']:
#         continue
#     genes.remove(gene_name.upper())
#     gene_synonyms[gene_name] = genes
#     # selected = unselected[(unselected['genesymbol_kinase'] == key) & (unselected['Kinases'].isin(genes))]
#     selected = unselected[(unselected['predictor'] == key) & (unselected['Kinases'].isin(genes))]
#     temp = pd.concat([temp, selected])
# endregion
#No additional overlaps were obtained.

#Filter out the predictions with context scores <0.9 and motif scores <0.5.
#There are 25 hits with large differences in the position numbers, so limit the difference to be 2
overlaps_NetworKIN['filter'] = overlaps_NetworKIN.apply(lambda x: (abs(x['Position_x']-x['Position_y'])<=2) and (x['context_score']>=0.9) and (x['motif_score']>=0.5), axis=1)
overlaps_NetworKIN = overlaps_NetworKIN[overlaps_NetworKIN['filter'] == True].iloc[: , :-1]
order = [0,3,2,1,4,9,10,7,8,11,5,6,12]
overlaps_NetworKIN = overlaps_NetworKIN[[overlaps_NetworKIN.columns[i] for i in order]]
overlaps_NetworKIN.drop_duplicates(inplace=True)
overlaps_NetworKIN.sort_values(by=['motif_score', 'context_score'], inplace=True, ascending=False)
# overlaps_NetworKIN.to_csv(r'..\Results\NetworKINvsPhospho.ELM.tsv', sep="\t", index=False)

#Perfect match
s3 = copy.deepcopy(overlaps_NetworKIN)
s3['filter'] = s3.apply(lambda x: x['Position_x']==x['Position_y'], axis=1)
s3 = s3[s3['filter'] == True].iloc[: , :-1]
order = [0,3,2,1,4,6,7,8,9,10,11,12]
s3 = s3[[s3.columns[i] for i in order]]
s3 = s3.rename(columns={s3.columns[1]: "Position"})
# s3.to_csv(r'..\Results\perfect_match_NetworKINvsPhospho.ELM.tsv', sep="\t", index=False)

#Confident overlaps
#Drop entries with missing "Kinases" and "KinaseID"
PhosphoELM_conf = PhosphoELM_format.dropna(subset=["KinaseID", "Kinases"], how="all")
#Map gene (family) names from "Kinases" to UniProt IDs
# region gname_to_uniprot
# genes= PhosphoELM_conf['Kinases'][PhosphoELM_conf['KinaseID'].isna()].unique().tolist()
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
# Add "Eg3 kinase" to gname_to_uniprot exclusively, as per the overlaps in "Phospho.ELM_vs_PhosphoSitePlus.py"
# Since "Eg3 kinase" is not a standard gene name, the format of query is without any field
# gene_name="Eg3 kinase"
# query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+{gene_name}&fields=accession'
# response = requests.get(query_url)
# uniprot = response.text.split('\n')[1:-1]
# gname_to_uniprot[gene_name] = uniprot
# with open('.\Variables\gname_to_uniprot.pkl', 'wb') as f:
#     pickle.dump(gname_to_uniprot, f)
PhosphoELM_conf = PhosphoELM_conf.copy() #To overcome "SettingWithCopyWarning" message
PhosphoELM_conf['KinaseID'] = PhosphoELM_conf.apply(lambda x: tuple(gname_to_uniprot[x['Kinases']]) if str(x['KinaseID'])=='nan' else x['KinaseID'], axis=1) #To overcome TypeError: unhashable type: 'list' while dropping duplicates
PhosphoELM_conf = PhosphoELM_conf[PhosphoELM_conf['KinaseID']!=()]
PhosphoELM_conf.drop('KinaseID_y', axis=1, inplace=True)
# Check for autophosphorylation
# PhosphoELM_conf['flag'] = NetworKIN_format.apply(lambda row: find_overlaps2(row['SubstrateID'], row['KinaseID']), axis=1)
# kk = PhosphoELM_conf.dropna(subset="flag") # No interactions
#Save as a pickle file to retain the data type of the tuple of strings
# with open('.\Variables\PhosphoELM_conf.pkl', 'wb') as f:
#     pickle.dump(PhosphoELM_conf, f)

#Drop the entries with "nan" in "Kinases" and Kinase IDs of Phospho.ELM
overlaps_NetworKIN_conf = overlaps_NetworKIN.dropna(subset=["KinaseID_y", "Kinases"], how="all")
overlaps_NetworKIN_conf = overlaps_NetworKIN_conf.loc[~(overlaps_NetworKIN_conf["Kinases"].str.contains("_group") & overlaps_NetworKIN_conf["KinaseID_y"].isna())]
overlaps_NetworKIN_conf.drop('KinaseID_y', axis=1, inplace=True)
overlaps_NetworKIN_conf = overlaps_NetworKIN_conf.rename(columns={overlaps_NetworKIN_conf.columns[1]: "KinaseID"})
overlaps_NetworKIN_conf.to_csv(r'..\Results\NetworKINvsPhospho.ELM.tsv', sep="\t", index=False)
s3_conf = s3.dropna(subset=["KinaseID_y", "Kinases"], how="all")
s3_conf = s3_conf.loc[~(s3_conf["Kinases"].str.contains("_group") & s3_conf["KinaseID_y"].isna())]
s3_conf.drop('KinaseID_y', axis=1, inplace=True)
s3_conf = s3_conf.rename(columns={s3_conf.columns[3]: "KinaseID"})
s3_conf.to_csv(r'..\Results\perfect_match_NetworKINvsPhospho.ELM.tsv', sep="\t", index=False)

#GPS and PhosphoSitePlus database
PhosphoSitePlus = pd.read_csv("..\Data\Kinase_Substrate_Dataset_formatted", sep="\t")
PhosphoSitePlus_subset = PhosphoSitePlus.loc[PhosphoSitePlus['KIN_ORGANISM'].eq('human') & PhosphoSitePlus['SUB_ORGANISM'].eq('human')] #Extract information for only homo sapiens
# PhosphoSitePlus_subset.isna().any()
#There are no missing values in the required columns
PhosphoSitePlus_format = PhosphoSitePlus_subset.iloc[:, [6,9,2,11]]
PhosphoSitePlus_format = PhosphoSitePlus_format.copy() #To overcome "SettingWithCopyWarning" message
PhosphoSitePlus_format['Code'] = PhosphoSitePlus_format["SUB_MOD_RSD"].apply(lambda x: str(x)[0])
PhosphoSitePlus_format['Position'] = PhosphoSitePlus_format["SUB_MOD_RSD"].apply(lambda x: int(str(x)[1:]))
PhosphoSitePlus_format.loc[:,'SITE_+/-7_AA']=PhosphoSitePlus_format['SITE_+/-7_AA'].apply(lambda x: x[:7].upper()+x[7]+x[8:].upper()) #Represent only the phosphorylation site in lowercase
# For positions less than 5, retain only the right end of the phosphorylation site (5-mer), and convert the rest to 9-mer
PhosphoSitePlus_format['SITE_+/-7_AA'] = PhosphoSitePlus_format.apply(lambda x: str(x['SITE_+/-7_AA'])[3:12] if x.Position > 5 else str(x['SITE_+/-7_AA'])[7:12], axis=1)
#If there is '_' in the 9-mer, slice it out
PhosphoSitePlus_format['SITE_+/-7_AA'] = PhosphoSitePlus_format['SITE_+/-7_AA'].str.replace('_', '', regex=True)
PhosphoSitePlus_format = PhosphoSitePlus_format.iloc[:, [0,5,4,2,3]]
PhosphoSitePlus_format.columns=['SubstrateID','Position','Code','KinaseID','Peptide']
# region kinase-kinase interactions
# #Extract only kinase-substrate interactions
# #Map uniprot ids of substrates to protein familes to filter out those which have "kinase"
# # region uniprot_to_pfam
# # uniprot= PhosphoSitePlus_format['SubstrateID'].unique().tolist()
# # base_url = 'https://rest.uniprot.org/uniprotkb'
# # uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# # for id in uniprot[2000:]:
# #     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
# #     response = requests.get(query_url, timeout=10)
# #     protein_fam = response.text.split('\n')[1]
# #     if "kinase" in protein_fam.lower():
# #         uniprot_to_pfam[id] = protein_fam
# # with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
# #     pickle.dump(uniprot_to_pfam, f)
# # endregion
# kinases = uniprot_to_pfam.keys()
# #Remove kinase-kinase interactions
# PhosphoSitePlus_format = PhosphoSitePlus_format[~PhosphoSitePlus_format['SubstrateID'].isin(kinases)]
# endregion
#Check for proteins that are not kinases in the kinase column
# region uniprot_to_pfam
# uniprot= PhosphoSitePlus_format['KinaseID'].unique().tolist()
# base_url = 'https://rest.uniprot.org/uniprotkb'
# uniprot = [x for x in uniprot if x not in uniprot_to_pfam.keys()]
# for id in uniprot:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     if "kinase" in protein_fam.lower():
#         uniprot_to_pfam[id] = protein_fam
# with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_pfam, f)
# endregion
# genes = PhosphoSitePlus[PhosphoSitePlus['KIN_ACC_ID'].isin(uniprot)]['GENE'].unique().tolist()
# genes = [x for x in genes if x not in gname_to_pfam.keys()]
#13 proteins that are not kinases are found - manually verified
#Exclude 'HSPA5' and 'ENPP3' and their corresponding UniProt IDs
# exclude = ['ENPP3', 'HSPA5']
# genes = [x for x in genes if x not in exclude]
# Update gname_to_pfam
# for gene_name in genes:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     gname_to_pfam[gene_name] = protein_fam
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# region uniprot_to_pfam
# ids = PhosphoSitePlus[PhosphoSitePlus['GENE'].isin(exclude)]['KIN_ACC_ID'].tolist()
# uniprot = [x for x in uniprot if x not in ids]
# for id in uniprot:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{id}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_fam = response.text.split('\n')[1]
#     uniprot_to_pfam[id] = protein_fam
# with open(r'.\Variables\uniprot_to_pfam.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_pfam, f)
# endregion
kinases = uniprot_to_pfam.keys()
PhosphoSitePlus_format = PhosphoSitePlus_format[PhosphoSitePlus_format['KinaseID'].isin(kinases)]
# Check for autophosphorylation
# kk = PhosphoSitePlus_format[PhosphoSitePlus_format['SubstrateID']==PhosphoSitePlus_format['KinaseID']] # No interactions
# PhosphoSitePlus_format.to_csv(r'.\Variables\PhosphoSitePlus_format.tsv', sep="\t", index=False)

overlaps_GPS = pd.merge(GPS_format, PhosphoSitePlus_format, how='inner', on=['SubstrateID','Code','Peptide'])
overlaps_GPS['KinaseID'] = overlaps_GPS.apply(lambda row: find_overlaps2(row['KinaseID_y'], row['KinaseID_x']), axis=1)
overlaps_GPS.dropna(subset="KinaseID", inplace=True)
overlaps_GPS.drop_duplicates(inplace=True)
#There are 81 hits with large differences in the position numbers, so limit the difference to be 2
overlaps_GPS['filter'] = overlaps_GPS.apply(lambda x: (abs(x['Position_x']-x['Position_y'])<=2), axis=1)
overlaps_GPS = overlaps_GPS[overlaps_GPS['filter'] == True].iloc[: , :-1]
order = [0,8,2,1,4,6,5]
overlaps_GPS = overlaps_GPS[[overlaps_GPS.columns[i] for i in order]]
overlaps_GPS.to_csv(r'..\Results\GPSvsPhosphoSitePlus.tsv', sep="\t", index=False)

#Perfect match
s4 = copy.deepcopy(overlaps_GPS)
s4['filter'] = s4.apply(lambda x: x['Position_x']==x['Position_y'], axis=1)
s4 = s4[s4['filter'] == True].iloc[: , :-1]
order = [0,3,2,1,4,6]
s4 = s4[[s4.columns[i] for i in order]]
s4 = s4.rename(columns={s4.columns[1]: "Position"})
s4.to_csv(r'..\Results\perfect_match_GPSvsPhosphoSitePlus.tsv', sep="\t", index=False)

##Number of vertices and edges
#vertices - kinases and substrates; edges - interaction between kinases and substrates (phosphorylation of substrates by kinases -> directed)
net_prop = pd.DataFrame({'Sources': ['GPS','NetworKIN','PhosphoSitePlus','Phospho.ELM'],
                         'Vertices': [len(GPS_format["KinaseID"].unique())+len(GPS_format["SubstrateID"].unique()),
                                      len(NetworKIN_format["KinaseID"].unique())+len(NetworKIN_format["SubstrateID"].unique()),
                                      len(PhosphoSitePlus_format['KinaseID'].unique())+len(PhosphoSitePlus_format['SubstrateID'].unique()),
                                      # len(PhosphoELM_format['Kinases'].unique())+len(PhosphoELM_format['SubstrateID'].unique())],
                                      len(PhosphoELM_conf['Kinases'].unique())+len(PhosphoELM_conf['SubstrateID'].unique())],
                         'Edges': [len(GPS_format.index), len(NetworKIN_format.index), len(PhosphoSitePlus_format.index), len(PhosphoELM_conf.index)]}).set_index("Sources")
# net_prop.to_csv(r'..\Results\network_properties.tsv', sep="\t")
#print(net_prop)

##Venn diagrams for the visualisations of the overlaps
#Labels are the number of edges (predictions)
matplotlib.use('Qt5Agg') #To prevent matplotlib from not responding
#Overlaps of predictions
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(GPS_format.index), len(NetworKIN_format.index), len(overlaps.index)),
             set_colors = ('purple', 'skyblue'),
             set_labels = ('GPS', 'NetworKIN'))
ax.set_title('GPS and NetworKIN \n (excluding phosphorylation site positions)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['GPS','Edges']-len(overlaps.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['NetworKIN','Edges']-len(overlaps.index)}")
venn.get_label_by_id('11').set_text(f"{len(overlaps.index)}")
venn.set_labels[1].set_x(0.35)
plt.legend([len(GPS_format.index), len(NetworKIN_format.index)],title="Total interactions",bbox_to_anchor=(0.8, 0.85))
plt.savefig('..\Results\Figures\overlaps_phosphorylation.png')

#Perfect match
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(GPS_format.index), len(NetworKIN_format.index), len(s2.index)),
             set_colors = ('purple', 'skyblue'),
             set_labels = ('GPS', 'NetworKIN'))
ax.set_title('GPS and NetworKIN \n (complete overlap)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['GPS','Edges']-len(s2.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['NetworKIN','Edges']-len(s2.index)}")
venn.get_label_by_id('11').set_text(f"{len(s2.index)}")
venn.set_labels[1].set_x(0.35)
plt.legend([len(GPS_format.index), len(NetworKIN_format.index)],title="Total interactions",bbox_to_anchor=(0.75, 0.85))
plt.savefig('..\Results\Figures\perfect_match_phosphorylation.png')

# #NetworKIN and Phospho.ELM database
# fig, ax = plt.subplots(figsize=(8, 6))
# venn = venn2(subsets = (len(PhosphoELM_format.index), len(NetworKIN_format.index), len(overlaps_NetworKIN.index)),
#              set_colors = ('orange', 'skyblue'),
#              set_labels = ('Phospho.ELM', 'NetworKIN'))
# ax.set_title('Phospho.ELM and NetworKIN \n (excluding phosphorylation site positions)')
# venn.get_patch_by_id('11').set_color('yellow')
# venn.get_label_by_id('10').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(overlaps_NetworKIN.index)}")
# venn.get_label_by_id('01').set_text(f"{net_prop.loc['NetworKIN','Edges']-len(overlaps_NetworKIN.index)}")
# venn.get_label_by_id('11').set_text(f"{len(overlaps_NetworKIN.index)}")
# plt.legend([len(PhosphoELM_format.index), len(NetworKIN_format.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.9))
# plt.savefig(r'..\Results\Figures\NetworKINvsPhospho.ELM.png') #Prefixed with 'r' to overcome the unicode error
#
# #NetworKIN and Phospho.ELM database - perfect match
# fig, ax = plt.subplots(figsize=(8, 6))
# venn = venn2(subsets = (len(PhosphoELM_format.index), len(NetworKIN_format.index), len(s3.index)),
#              set_colors = ('orange', 'skyblue'),
#              set_labels = ('Phospho.ELM', 'NetworKIN'))
# ax.set_title('Phospho.ELM and NetworKIN \n (complete overlap)')
# venn.get_patch_by_id('11').set_color('yellow')
# venn.get_label_by_id('10').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(s3.index)}")
# venn.get_label_by_id('01').set_text(f"{net_prop.loc['NetworKIN','Edges']-len(s3.index)}")
# venn.get_label_by_id('11').set_text(f"{len(s3.index)}")
# plt.legend([len(PhosphoELM_format.index), len(NetworKIN_format.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.9))
# plt.savefig('..\Results\Figures\perfect_match_NetworKINvsPhospho.ELM.png')

#NetworKIN and Phospho.ELM database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PhosphoELM_conf.index), len(NetworKIN_format.index), len(overlaps_NetworKIN_conf.index)),
             set_colors = ('orange', 'skyblue'),
             set_labels = ('Phospho.ELM', 'NetworKIN'))
ax.set_title('Phospho.ELM and NetworKIN \n (excluding phosphorylation site positions)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(overlaps_NetworKIN_conf.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['NetworKIN','Edges']-len(overlaps_NetworKIN_conf.index)}")
venn.get_label_by_id('11').set_text(f"{len(overlaps_NetworKIN_conf.index)}")
venn.set_labels[1].set_x(0.03)
venn.set_labels[0].set_x(-0.4)
plt.legend([len(PhosphoELM_conf.index), len(NetworKIN_format.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.91))
plt.savefig(r'..\Results\Figures\NetworKINvsPhospho.ELM.png') #Prefixed with 'r' to overcome the unicode error

#NetworKIN and Phospho.ELM database - perfect match
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PhosphoELM_conf.index), len(NetworKIN_format.index), len(s3_conf.index)),
             set_colors = ('orange', 'skyblue'),
             set_labels = ('Phospho.ELM', 'NetworKIN'))
ax.set_title('Phospho.ELM and NetworKIN \n (complete overlap)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(s3_conf.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['NetworKIN','Edges']-len(s3_conf.index)}")
venn.get_label_by_id('11').set_text(f"{len(s3_conf.index)}")
venn.set_labels[1].set_x(0.03)
venn.set_labels[0].set_x(-0.4)
plt.legend([len(PhosphoELM_conf.index), len(NetworKIN_format.index)],title="Total interactions",bbox_to_anchor=(0.75, 0.91))
plt.savefig('..\Results\Figures\perfect_match_NetworKINvsPhospho.ELM.png')

#GPS and PhosphoSitePlus database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PhosphoSitePlus_format.index), len(GPS_format.index), len(overlaps_GPS.index)),
             set_colors = ('orange', 'purple'),
             set_labels = ('PhosphoSitePlus', 'GPS'))
ax.set_title('PhosphoSitePlus and GPS \n (excluding phosphorylation site positions)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(overlaps_GPS.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['GPS','Edges']-len(overlaps_GPS.index)}")
venn.get_label_by_id('11').set_text(f"{len(overlaps_GPS.index)}")
plt.legend([len(PhosphoSitePlus_format.index), len(GPS_format.index)],title="Total interactions",bbox_to_anchor=(0.85, 0.9))
plt.savefig('..\Results\Figures\GPSvsPhosphoSitePlus.png')

#GPS and PhosphoSitePlus database - perfect match
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PhosphoSitePlus_format.index), len(GPS_format.index), len(overlaps_GPS.index)),
             set_colors = ('orange', 'purple'),
             set_labels = ('PhosphoSitePlus', 'GPS'))
ax.set_title('PhosphoSitePlus and GPS \n (complete overlap)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(s4.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['GPS','Edges']-len(s4.index)}")
venn.get_label_by_id('11').set_text(f"{len(s4.index)}")
plt.legend([len(PhosphoSitePlus_format.index), len(GPS_format.index)],title="Total interactions",bbox_to_anchor=(0.75, 0.9))
plt.savefig('..\Results\Figures\perfect_match_GPSvsPhosphoSitePlus.png')

##Sensitivity
#NetworKIN and Phospho.ELM
true_positives = len(overlaps_NetworKIN_conf)
false_negatives = len(PhosphoELM_conf) - true_positives
sensitivity = true_positives / (true_positives + false_negatives)
print("Sensitivity (NetworKIN):", sensitivity)

#GPS and PhosphoSitePlus
true_positives = len(overlaps_GPS)
false_negatives = len(PhosphoSitePlus_format) - true_positives
sensitivity = true_positives / (true_positives + false_negatives)
print("Sensitivity (GPS):", sensitivity)