import pandas as pd
import copy
import pickle
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import matplotlib
import requests

#Load saved variables
PhosphoELM_format = pd.read_csv(r'.\Variables\PhosphoELM_format.tsv', sep="\t")
with open('.\Variables\gname_to_uniprot.pkl', 'rb') as f:
    gname_to_uniprot = pickle.load(f)
with open('.\Variables\kin_dict.pkl', 'rb') as f:
    kin_dict = pickle.load(f)
with open(r'.\Variables\not_in_kin_dict.pkl', 'rb') as f:
    not_in_kin_dict = pickle.load(f)
with open(r'.\Variables\uniprot_to_pfam.pkl', 'rb') as f:
    uniprot_to_pfam = pickle.load(f)
net_prop = pd.read_csv(r'..\Results\network_properties.tsv', sep="\t", index_col=0)

#Format PhosphoSitePlus with the addition of the column "KINASE"
PhosphoSitePlus = pd.read_csv("..\Data\Kinase_Substrate_Dataset_formatted", sep="\t")
PhosphoSitePlus_subset = PhosphoSitePlus.loc[PhosphoSitePlus['KIN_ORGANISM'].eq('human') & PhosphoSitePlus['SUB_ORGANISM'].eq('human')] #Extract information for only homo sapiens
# PhosphoSitePlus_subset.isna().any()
#There are no missing values in the required columns
PhosphoSitePlus_format = PhosphoSitePlus_subset.iloc[:, [6,9,2,11,0,1]]
PhosphoSitePlus_format = PhosphoSitePlus_format.copy() #To overcome "SettingWithCopyWarning" message
PhosphoSitePlus_format['Code'] = PhosphoSitePlus_format["SUB_MOD_RSD"].apply(lambda x: str(x)[0])
PhosphoSitePlus_format['Position'] = PhosphoSitePlus_format["SUB_MOD_RSD"].apply(lambda x: int(str(x)[1:]))
PhosphoSitePlus_format.loc[:,'SITE_+/-7_AA']=PhosphoSitePlus_format['SITE_+/-7_AA'].apply(lambda x: x[:7].upper()+x[7]+x[8:].upper()) #Represent only the phosphorylation site in lowercase
# For positions less than 5, retain only the right end of the phosphorylation site (5-mer), and convert the rest to 9-mer
PhosphoSitePlus_format['SITE_+/-7_AA'] = PhosphoSitePlus_format.apply(lambda x: str(x['SITE_+/-7_AA'])[3:12] if x.Position > 5 else str(x['SITE_+/-7_AA'])[7:12], axis=1)
#If there is '_' in the 9-mer, slice it out
PhosphoSitePlus_format['SITE_+/-7_AA'] = PhosphoSitePlus_format['SITE_+/-7_AA'].str.replace('_', '', regex=True)
PhosphoSitePlus_format = PhosphoSitePlus_format.iloc[:, [0,7,6,2,3,5,4]]
PhosphoSitePlus_format.columns=['SubstrateID','Position','Code','KinaseID','Peptide']+PhosphoSitePlus_format.columns[5:].tolist()
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
kinases = uniprot_to_pfam.keys()
PhosphoSitePlus_format = PhosphoSitePlus_format[PhosphoSitePlus_format['KinaseID'].isin(kinases)]

##Overlaps
#1. Overlaps for the entries having a UniProt ID for kinases
selected = PhosphoELM_format[~PhosphoELM_format['KinaseID'].isna()]
unselected = PhosphoELM_format[PhosphoELM_format['KinaseID'].isna()]
overlaps = pd.merge(PhosphoSitePlus_format, selected, how='inner', on=['SubstrateID','KinaseID','Code','Peptide'])
#2. Overlaps for the entries without UniProt IDs for kinases
unselected= unselected.copy() #To overcome "SettingWithCopyWarning" message
unselected.drop('KinaseID_y',axis=1, inplace=True)
overlaps_unselected = pd.merge(PhosphoSitePlus_format, unselected, how='inner', on=['SubstrateID','Code','Peptide'])
unselected = overlaps_unselected[~overlaps_unselected['Kinases'].isna()]
selected = overlaps_unselected[overlaps_unselected['Kinases'].isna()]
overlaps.columns=unselected.columns #Rename "KinaseID" as "KinaseID_x"
overlaps = pd.concat([overlaps, selected])
#3. Manually inspect the 'GENE','KINASE' and 'Kinases' columns of the unselected entries.
#'GENE'
#Select entries where "GENE" and "Kinases" are the same (only 6 such entries)
selected = unselected.query('GENE==Kinases')
overlaps = pd.concat([overlaps, selected])
unselected = unselected.query('GENE!=Kinases')
#Use previous curation of 'kin_dict' first
for key in kin_dict:
    kinase_list = kin_dict.get(key)
    selected = unselected[(unselected['GENE'] == key) & (unselected['Kinases'].isin(kinase_list))]
    overlaps = pd.concat([overlaps, selected])
mask = ~unselected[unselected.columns].apply(tuple, axis=1).isin(overlaps[overlaps.columns].apply(tuple, axis=1))
unselected = unselected[mask]
#New information
all_kin_dict = {}
for predictor in unselected["GENE"].unique():
    kinases = unselected.loc[unselected["GENE"] == predictor, "Kinases"].unique()
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
new_info = {"PRKCD": ["PKC_group"],
            "MAPK8": ["MAPK_group","JNK_group"],
            "PRKCB": ["PKC_group"],
            "PRKACA": ["PKA_group"],
            "PRKCH": ["PKC_group"],
            "PRKCZ": ["PKC_group"],
            "AKT3": ["PKB_group"],
            "RPS6KA1": ["RSK_group"],
            "RPS6KA3": ["RSK_group"],
            "CSNK2A1": ["CK2_group"],
            "MAPK3": ["MAPK_group"],
            "MAPK14": ["p38_group"],
            "PRKD1": ["PKG/cGK_group"],
            "PRKCE": ["PKC_group"],
            "ROCK2": ["ROCK_group"],
            "PRKAA1": ["AMPK_group"],
            "PKN1": ["PKC_group"],
            "PKN3": ["PKC_group"],
            "IKBKE": ["IKK_group"],
            "SGK1": ["SGK_group"],
            "RPS6KB2": ["p70S6K_group", "RSK_group"],
            "PRKCI": ["PKC_group"],
            "CDK1": ["CDK_group"],
            "CAMK2A": ["CaM-KII_group"],
            "CDK6": ["CDK_group"],
            "RPS6KA5": ["RSK_group"],
            "GRK2": ["GRK_group"],
            "DAPK1": ["DAPK_group"],
            "AKT2": ["PKB_group"],
            "MARK2": ["MARK_group"],
            "GSK3A": ["GSK-3_group"],
            "CDK5": ["CDK_group"],
            "CSNK1A1": ["CK1_group"],
            "MAP2K7": ["MAP2K_group"],
            "CHUK": ["IKK_group"],
            "IKBKB": ["IKK_group"],
            "NLK": ["MAPK_group"],
            "GRK5": ["GRK_group"],
            "ROCK1": ["ROCK_group"],
            "PAK1": ["PAK_group"],
            "PDGFRB": ["PDGFR_group"],
            "CDK7": ["CDK_group"],
            "CDK3": ["CDK_group"],
            "MAPK12": ["MAPK_group"],
            "CDK9": ["CDK_group"],
            "CDK18": ["CDK_group"],
            "CDK4": ["CDK_group"],
            "PRKAA2": ["AMPK_group"],
            "RPS6KA4": ["RSK_group"],
            "MAP3K8": ["MAP3K_group"],
            "RPS6KA2": ["RSK_group"],
            "TBK1": ["IKK_group"],
            "CSNK1E": ["CK1_group"],
            "SRC": ["SRC_group"],
            "FYN": ["SRC_group"],
            "MAPK11": ["MAPK_group"],
            "MAPK13": ["MAPK_group"],
            "MAP2K1": ["MAP2K_group"],
            "PAK2": ["PAK_group"],
            "MELK": ["Eg3 kinase"],
            "PRKACB": ["PKA_group"],
            "DMPK": ["DMPK_group"],
            "MAPK10": ["MAPK_group"],
            "MAP2K2": ["MAP2K_group"],
            "FGFR3": ["FGFR_group"],
            "GRK6": ["GRK_group"],
            "GRK3": ["GRK_group"],
            "MAP2K4": ["MAP2K_group"],
            "MAP2K6": ["MAP2K_group"],
            "MAP2K3": ["MAP2K_group"],
            "LCK": ["SRC_group"],
            "HCK": ["SRC_group"],
            "MAPK7": ["MAPK_group"],
            "MAPK6": ["MAPK_group"],
            "PRKG1": ["PKG/cGK_group"],
            "PHKA1": ["PHK_group"]}
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
    selected = unselected[(unselected['GENE'] == key) & (unselected['Kinases'].isin(kinase_list))]
    overlaps = pd.concat([overlaps, selected])
mask = ~unselected[unselected.columns].apply(tuple, axis=1).isin(overlaps[overlaps.columns].apply(tuple, axis=1))
unselected = unselected[mask]
#'KINASE'
#Use previous curation of 'kin_dict' first
for key in kin_dict:
    kinase_list = kin_dict.get(key)
    selected = unselected[(unselected['KINASE'] == key) & (unselected['Kinases'].isin(kinase_list))]
    overlaps = pd.concat([overlaps, selected])
mask = ~unselected[unselected.columns].apply(tuple, axis=1).isin(overlaps[overlaps.columns].apply(tuple, axis=1))
unselected = unselected[mask]
#New information
all_kin_dict = {}
for predictor in unselected["KINASE"].unique():
    kinases = unselected.loc[unselected["KINASE"] == predictor, "Kinases"].unique()
    all_kin_dict[predictor] = list(kinases)
# df = pd.DataFrame({'key' : all_kin_dict.keys() , 'value' : all_kin_dict.values() })
# df.to_csv(r'.\Variables\all_kin_dict.tsv', sep="\t", index=False)
#No new information
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
# with open('.\Variables\kin_dict.pkl', 'wb') as f:
#     pickle.dump(kin_dict, f)
# with open(r'.\Variables\not_in_kin_dict.pkl', 'wb') as f:
#     pickle.dump(not_in_kin_dict, f)

# region Synonyms
# #Get the synonyms of kinases from UniProtKB (Gene Names column)
# base_url = 'https://rest.uniprot.org/uniprotkb'
# # gene_names = unselected["GENE"].unique()
# gene_names = unselected["KINASE"].unique()
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
#     # selected = unselected[(unselected['GENE'] == key) & (unselected['Kinases'].isin(genes))]
#     selected = unselected[(unselected['KINASE'] == key) & (unselected['Kinases'].isin(genes))]
#     temp = pd.concat([temp, selected])
# endregion
#No additional overlaps were obtained.

#Set the difference in position number to a maximum of 3. Beyond that, most of the "Kinases" have 'nan' or the position numbers have large differences.
overlaps['filter'] = overlaps.apply(lambda x: abs(x['Position_x']-x['Position_y'])<=3, axis=1)
overlaps = overlaps[overlaps['filter'] == True].iloc[: , :-1]
order = [0,3,2,1,4,7,8,5,6,9,10]
overlaps = overlaps[[overlaps.columns[i] for i in order]]
overlaps.drop_duplicates(inplace=True)
# overlaps.to_csv(r'..\Results\Phospho.ELMvsPhosphoSitePlus.tsv', sep="\t", index=False)

#Perfect match
s1 = copy.deepcopy(overlaps)
s1['filter'] = s1.apply(lambda x: x['Position_x']==x['Position_y'], axis=1)
s1 = s1[s1['filter'] == True].iloc[: , :-1]
order = [0,3,2,1,4,6,7,8,9,10]
s1 = s1[[s1.columns[i] for i in order]]
s1 = s1.rename(columns={s1.columns[1]: "Position"})
# s1.to_csv(r'..\Results\perfect_match_Phospho.ELMvsPhosphoSitePlus.tsv', sep="\t", index=False)

#Confident overlaps
PhosphoELM_format["KinaseID"] = PhosphoELM_format["KinaseID"].apply(lambda x: tuple(x.split(",")) if str(x)!="nan" else x)
#Drop entries with missing "Kinases" and "KinaseID"
PhosphoELM_conf = PhosphoELM_format.dropna(subset=["KinaseID", "Kinases"], how="all")
#Map gene (family) names from "Kinases" to UniProt IDs
# genes= PhosphoELM_conf['Kinases'][PhosphoELM_conf['KinaseID'].isna()].unique().tolist()
# gene_names = [x for x in genes if x not in gname_to_uniprot.keys()]
#No additional gene names to be added to "gname_to_uniprot"
PhosphoELM_conf = PhosphoELM_conf.copy() #To overcome "SettingWithCopyWarning" message
PhosphoELM_conf['KinaseID'] = PhosphoELM_conf.apply(lambda x: tuple(gname_to_uniprot[x['Kinases']]) if str(x['KinaseID'])=='nan' else x['KinaseID'], axis=1) #To overcome TypeError: unhashable type: 'list' while dropping duplicates
PhosphoELM_conf = PhosphoELM_conf[PhosphoELM_conf['KinaseID']!=()]
PhosphoELM_conf.drop('KinaseID_y', axis=1, inplace=True)
# with open('.\Variables\PhosphoELM_conf.pkl', 'wb') as f:
#     pickle.dump(PhosphoELM_conf, f)

#Drop the entries with "nan" in "Kinases" and Kinase IDs of Phospho.ELM
overlaps_conf = overlaps.dropna(subset=["KinaseID_y", "Kinases"], how="all")
overlaps_conf = overlaps_conf.loc[~(overlaps_conf["Kinases"].str.contains("_group") & overlaps_conf["KinaseID_y"].isna())]
overlaps_conf.drop('KinaseID_y', axis=1, inplace=True)
overlaps_conf = overlaps_conf.rename(columns={overlaps_conf.columns[1]: "KinaseID"})
overlaps_conf.to_csv(r'..\Results\Phospho.ELMvsPhosphoSitePlus.tsv', sep="\t", index=False)
s1_conf = s1.dropna(subset=["KinaseID_y", "Kinases"], how="all")
s1_conf = s1_conf.loc[~(s1_conf["Kinases"].str.contains("_group") & s1_conf["KinaseID_y"].isna())]
s1_conf.drop('KinaseID_y', axis=1, inplace=True)
s1_conf = s1_conf.rename(columns={s1_conf.columns[3]: "KinaseID"})
s1_conf.to_csv(r'..\Results\perfect_match_Phospho.ELMvsPhosphoSitePlus.tsv', sep="\t", index=False)

##Venn diagrams for the visualisations of the overlaps
#Labels are the number of edges (predictions)
matplotlib.use('Qt5Agg') #To prevent matplotlib from not responding
# #PhosphoSitePlus and Phospho.ELM database
# fig, ax = plt.subplots(figsize=(8, 6))
# venn = venn2(subsets = (len(PhosphoELM_format.index), len(PhosphoSitePlus_format.index), len(overlaps.index)),
#              set_colors = ('brown', 'orange'),
#              set_labels = ('Phospho.ELM', 'PhosphoSitePlus'))
# ax.set_title('Phospho.ELM and PhosphoSitePlus \n (excluding phosphorylation site positions)')
# venn.get_patch_by_id('11').set_color('yellow')
# venn.get_label_by_id('01').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(overlaps.index)}")
# venn.get_label_by_id('10').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(overlaps.index)}")
# venn.get_label_by_id('11').set_text(f"{len(overlaps.index)}")
# plt.legend([len(PhosphoELM_format.index), len(PhosphoSitePlus_format.index)],title="Total interactions")
# plt.savefig(r'..\Results\Figures\Phospho.ELMvsPhosphoSitePlus.png') #Prefixed with 'r' to overcome the unicode error
#
# #PhosphoSitePlus and Phospho.ELM database - perfect match
# fig, ax = plt.subplots(figsize=(8, 6))
# venn = venn2(subsets = (len(PhosphoELM_format.index), len(PhosphoSitePlus_format.index), len(s1.index)),
#              set_colors = ('brown','orange'),
#              set_labels = ('Phospho.ELM', 'PhosphoSitePlus'))
# ax.set_title('Phospho.ELM and PhosphoSitePlus \n (complete overlap)')
# venn.get_patch_by_id('11').set_color('yellow')
# venn.get_label_by_id('01').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(s1.index)}")
# venn.get_label_by_id('10').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(s1.index)}")
# venn.get_label_by_id('11').set_text(f"{len(s1.index)}")
# plt.legend([len(PhosphoELM_format.index), len(PhosphoSitePlus_format.index)],title="Total interactions")
# plt.savefig(r'..\Results\Figures\perfect_match_Phospho.ELMvsPhosphoSitePlus.png')

#PhosphoSitePlus and Phospho.ELM database
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PhosphoSitePlus_format.index), len(PhosphoELM_conf.index), len(overlaps_conf.index)),
             set_colors = ('orange','brown'),
             set_labels = ('PhosphoSitePlus','Phospho.ELM'))
ax.set_title('Phospho.ELM and PhosphoSitePlus \n (excluding phosphorylation site positions)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(overlaps_conf.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(overlaps_conf.index)}")
venn.get_label_by_id('11').set_text(f"{len(overlaps_conf.index)}")
plt.legend([len(PhosphoSitePlus_format.index),len(PhosphoELM_conf.index)],title="Total interactions")
plt.savefig(r'..\Results\Figures\Phospho.ELMvsPhosphoSitePlus.png') #Prefixed with 'r' to overcome the unicode error

#PhosphoSitePlus and Phospho.ELM database - perfect match
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PhosphoSitePlus_format.index), len(PhosphoELM_conf.index), len(s1_conf.index)),
             set_colors = ('orange','brown'),
             set_labels = ('PhosphoSitePlus','Phospho.ELM'))
ax.set_title('Phospho.ELM and PhosphoSitePlus \n (complete overlap)')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(s1_conf.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(s1_conf.index)}")
venn.get_label_by_id('11').set_text(f"{len(s1_conf.index)}")
plt.legend([len(PhosphoSitePlus_format.index),len(PhosphoELM_conf.index)],title="Total interactions")
plt.savefig(r'..\Results\Figures\perfect_match_Phospho.ELMvsPhosphoSitePlus.png')