import pandas as pd
import pickle
from modules import *
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Patch
import requests

#Load saved variables
net_prop = pd.read_csv(r'..\Results\network_properties.tsv', sep="\t", index_col=0)
with open(r'.\Variables\uniprot_to_gname.pkl', 'rb') as f:
    uniprot_to_gname = pickle.load(f)

##Read the input files of the databases and interactomes
with open('.\Variables\PhosphoELM_conf.pkl', 'rb') as f:
    PhosphoELM_conf = pickle.load(f)
PhosphoSitePlus_format = pd.read_csv(r'.\Variables\PhosphoSitePlus_format.tsv', sep="\t")
STRING_format = pd.read_csv(r'.\Variables\STRING_format.tsv', sep="\t")
with open('.\Variables\BioGRID_format.pkl', 'rb') as f:
    BioGRID_format = pickle.load(f)
with open('.\Variables\PC_format.pkl', 'rb') as f:
    PC_format = pickle.load(f)
PhosphoELMvsPhosphoSitePlus_conf = pd.read_csv(r'..\Results\Phospho.ELMvsPhosphoSitePlus.tsv', sep="\t")

##Phospho.ELM vs Interactomes
#STRING
merge1 = check_condition2(STRING_format, PhosphoELM_conf)
merge1 = merge1.reindex(columns=list(PhosphoELM_conf.columns) + list(STRING_format.columns))
merge1['KinaseID'] = merge1['protein2']
merge1.drop_duplicates(inplace=True)
merge2 = check_condition2_reverse(STRING_format, PhosphoELM_conf)
merge2 = merge2.reindex(columns=list(PhosphoELM_conf.columns) + list(STRING_format.columns))
merge2['KinaseID'] = merge2['protein1']
merge2.drop_duplicates(inplace=True)
PhosphoELM_STRING = pd.concat([merge1, merge2])
PhosphoELM_STRING = PhosphoELM_STRING.drop(['protein1', 'protein2'], axis=1)
PhosphoELM_STRING.to_csv(r'..\Results\Phospho.ELMvsSTRING.tsv', sep="\t", index=False)

#BioGRID
#Extract rows that have a single element in the "protein1" column and convert them to string
single_pr1 = BioGRID_format[BioGRID_format['protein1'].apply(lambda x: len(x) == 1)]
single_pr1 = single_pr1.copy()
single_pr1['protein1'] = single_pr1['protein1'].apply(lambda x: str(x[0]))
# Extract rows with multiple elements in the "protein1" tuple
mul_pr1 = BioGRID_format[BioGRID_format['protein1'].apply(lambda x: len(x) > 1)]
merge1_single = check_condition3(single_pr1, PhosphoELM_conf)
merge1_single = merge1_single.reindex(columns=list(PhosphoELM_conf.columns) + list(single_pr1.columns))
merge1_mul = check_condition4(mul_pr1, PhosphoELM_conf)
#No interactions - verified manually

#Extract rows that have a single element in the "protein2" column and convert them to string
single_pr2 = BioGRID_format[BioGRID_format['protein2'].apply(lambda x: len(x) == 1)]
single_pr2 = single_pr2.copy()
single_pr2['protein2'] = single_pr2['protein2'].apply(lambda x: str(x[0]))
# Extract rows with multiple elements in the "protein2" tuple
mul_pr2 = BioGRID_format[BioGRID_format['protein2'].apply(lambda x: len(x) > 1)]
merge2_single = check_condition3_reverse(single_pr2, PhosphoELM_conf)
merge2_single = merge2_single.reindex(columns=list(PhosphoELM_conf.columns) + list(single_pr2.columns))
merge2_mul = check_condition4_reverse(mul_pr2, PhosphoELM_conf)
#No interactions - verified manually

PhosphoELM_BioGRID = pd.concat([merge1_single, merge2_single])
PhosphoELM_BioGRID = PhosphoELM_BioGRID.drop(['protein1', 'protein2'], axis=1)
with open('..\Results\Phospho.ELMvsBioGRID.pkl', 'wb') as f:
    pickle.dump(PhosphoELM_BioGRID, f)

#Pathway Commons
merge1 = check_condition4(PC_format, PhosphoELM_conf)
merge1 = merge1.reindex(columns=list(PhosphoELM_conf.columns) + list(PC_format.columns))
merge2 = check_condition4_reverse(PC_format, PhosphoELM_conf)
merge2 = merge2.reindex(columns=list(PhosphoELM_conf.columns) + list(PC_format.columns))
PhosphoELM_PC = pd.concat([merge1, merge2])
PhosphoELM_PC = PhosphoELM_PC.drop(['protein1', 'protein2'], axis=1)
PhosphoELM_PC.drop_duplicates(inplace=True)
with open('..\Results\Phospho.ELMvsPC.pkl', 'wb') as f:
    pickle.dump(PhosphoELM_PC, f)

##PhosphoSitePlus vs Interactomes
#STRING
merge1 = pd.merge(PhosphoSitePlus_format, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein1', 'protein2'])
merge2 = pd.merge(PhosphoSitePlus_format, STRING_format, left_on=['SubstrateID', 'KinaseID'], right_on=['protein2', 'protein1'])
PhosphoSitePlus_STRING = pd.concat([merge1, merge2])
PhosphoSitePlus_STRING = PhosphoSitePlus_STRING.drop(['protein1', 'protein2'], axis=1)
PhosphoSitePlus_STRING.to_csv(r'..\Results\PhosphoSitePlusvsSTRING.tsv', sep="\t", index=False)

#BioGRID
merge1 = pd.concat(BioGRID_format.apply(check_condition, args=(PhosphoSitePlus_format,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(PhosphoSitePlus_format.columns) + list(BioGRID_format.columns))
merge2 = pd.concat(BioGRID_format.apply(check_condition_reverse, args=(PhosphoSitePlus_format,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(PhosphoSitePlus_format.columns) + list(BioGRID_format.columns))
PhosphoSitePlus_BioGRID = pd.concat([merge1, merge2])
PhosphoSitePlus_BioGRID = PhosphoSitePlus_BioGRID.drop(['protein1', 'protein2'], axis=1)
PhosphoSitePlus_BioGRID.to_csv(r'..\Results\PhosphoSitePlusvsBioGRID.tsv', sep="\t", index=False)

#Pathway Commons
merge1 = pd.concat(PC_format.apply(check_condition, args=(PhosphoSitePlus_format,), axis=1).tolist(), ignore_index=True)
merge1 = merge1.reindex(columns=list(PhosphoSitePlus_format.columns) + list(PC_format.columns))
merge2 = pd.concat(PC_format.apply(check_condition_reverse, args=(PhosphoSitePlus_format,), axis=1).tolist(), ignore_index=True)
merge2 = merge2.reindex(columns=list(PhosphoSitePlus_format.columns) + list(PC_format.columns))
PhosphoSitePlus_PC = pd.concat([merge1, merge2])
PhosphoSitePlus_PC = PhosphoSitePlus_PC.drop(['protein1', 'protein2'], axis=1)
#No duplicates
PhosphoSitePlus_PC.to_csv(r'..\Results\PhosphoSitePlusvsPC.tsv', sep="\t", index=False)

# region BioGRID vs BioGRID_PC
PC_BioGRID = PC_format[PC_format['INTERACTION_DATA_SOURCE'].str.contains("BioGRID")]
merge1 = pd.merge(BioGRID_format, PC_BioGRID, left_on=['Official Symbol Interactor A', 'Official Symbol Interactor B'], right_on=['PARTICIPANT_A', 'PARTICIPANT_B'])
merge2 = pd.merge(BioGRID_format, PC_BioGRID, left_on=['Official Symbol Interactor A', 'Official Symbol Interactor B'], right_on=['PARTICIPANT_B', 'PARTICIPANT_A'])
overlaps = pd.concat([merge1, merge2])
# There are 2717 interactions in the BioGRID database of PC that are not in the BioGRID database.
# There are 51166 interactions in the BioGRID database that are not classified in PC.

# Check if nodes of the exclusive interactions of the 2019 version is present in the latest version
overlaps.rename(columns={'protein1_y': 'protein1', 'protein2_y': 'protein2'}, inplace=True)
overlaps1 = overlaps.iloc[:,[4,5,6,7,8,9]]
overlaps1=overlaps1.copy()
PC_BioGRID_excl = pd.merge(PC_BioGRID, overlaps1, how='left', indicator=True)
PC_BioGRID_excl = PC_BioGRID_excl[PC_BioGRID_excl['_merge'] == 'left_only'].drop(columns=['_merge'])
node_PC = pd.concat([PC_BioGRID['PARTICIPANT_A'], PC_BioGRID['PARTICIPANT_B']]).unique()
node_BioGRID = pd.concat([BioGRID_format['Official Symbol Interactor A'], BioGRID_format['Official Symbol Interactor B']]).unique()
excl_nodes = [x for x in node_PC if x not in node_BioGRID]
# There are 341 exclusive nodes in the BioGRID 2019 version of PC that are not in the latest BioGRID version

# Rate of increase in publications per month
# print(((83750-83418)+(83418-83191)+(83191-82951)+(82951-82750)+(82750-82612)+(82024-81882)+(82134-82024)+(82354-82134)+(82612-82354))/10) #187
# endregion

##Classify Pathway Commons to its individual databases
# Position_x - PhosphoSitePlus, Position_y - Phospho.ELM
# Interactions that are exclusively present in Phospho.ELM
merged = pd.merge(PhosphoELM_conf, PhosphoELMvsPhosphoSitePlus_conf, how='inner', on=['SubstrateID', 'Code', 'Peptide'])
merged  = merged[merged ['Position'] == merged['Position_y']]
merged['KinaseID'] = merged.apply(lambda row: find_overlaps2(row['KinaseID_y'], row['KinaseID_x']),axis=1)
merged.dropna(subset="KinaseID", inplace=True)
merged = merged.iloc[:, :7]
merged.rename(columns={'KinaseID_x': 'KinaseID', 'Kinases_x': 'Kinases', 'Description_x': 'Description'}, inplace=True)
PhosphoELM_exl = PhosphoELM_conf.merge(merged, how='left', indicator=True)
PhosphoELM_exl = PhosphoELM_exl[PhosphoELM_exl['_merge'] == 'left_only']
PhosphoELM_exl = PhosphoELM_exl.drop(columns=['_merge'])
# Interactions that are exclusively present in PhosphoSitePlus
merged = pd.merge(PhosphoSitePlus_format, PhosphoELMvsPhosphoSitePlus_conf, how='inner', on=['SubstrateID', 'Code', 'Peptide', 'KinaseID'])
merged  = merged[merged ['Position'] == merged['Position_x']]
merged = merged.iloc[:, :5]
PhosphoSitePlus_exl = PhosphoSitePlus_format.merge(merged, how='left', indicator=True)
PhosphoSitePlus_exl = PhosphoSitePlus_exl[PhosphoSitePlus_exl['_merge'] == 'left_only']
PhosphoSitePlus_exl = PhosphoSitePlus_exl.drop(columns=['_merge'])

PC_databases = PC_format['INTERACTION_DATA_SOURCE'].str.split(';').explode('0').unique().tolist()
PC_prop=pd.DataFrame()
for i in PC_databases:
    PhosphoELM_i = PhosphoELM_PC[PhosphoELM_PC['INTERACTION_DATA_SOURCE'].str.contains(i)]
    PhosphoSitePlus_i = PhosphoSitePlus_PC[PhosphoSitePlus_PC['INTERACTION_DATA_SOURCE'].str.contains(i)]
    #Retain only unique interactions
    PhosphoELM_i = PhosphoELM_i.drop_duplicates(subset=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'])
    PhosphoSitePlus_i = PhosphoSitePlus_i.drop_duplicates(subset=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'])
    merge1 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoELM_i, how='inner', on=['SubstrateID', 'Code', 'Peptide'])
    if not merge1.empty:
        merge1['KinaseID'] = merge1.apply(lambda row: find_overlaps2(row['KinaseID_x'], row['KinaseID_y']),axis=1)
        merge1.dropna(subset="KinaseID", inplace=True)
        merge1 = merge1.iloc[:, [0, 18, 2, 3, 4, 5]]
    else:
        merge1 = pd.DataFrame()
    merge2 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoSitePlus_i, how='inner', on=['SubstrateID','Code','KinaseID','Peptide'])
    merge2 = merge2.iloc[:, [0, 1, 2, 3, 4, 5]]
    overlaps_conf = pd.concat([merge1, merge2])
    overlaps_conf.drop_duplicates(inplace=True)
    # High confidence additions
    merged = PhosphoELMvsPhosphoSitePlus_conf.merge(overlaps_conf, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide','Position_x','Position_y'], how='left', indicator=True)
    # Filter rows that are only in Phospho.ELM and PhosphoSitePlus, but not in interactome
    high_conf = merged[merged['_merge'] == 'left_only']
    high_conf = high_conf.drop(columns=['_merge'])
    high_conf['Confidence'] = 'High'
    high_conf = high_conf.iloc[:, :6].join(high_conf.iloc[:, -1])
    # Medium confidence additions
    # Phospho.ELM
    merged = PhosphoELM_exl.merge(PhosphoELM_i, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'], how='left', indicator=True)
    # Filter rows that are only in Phospho.ELM, but not in interactome
    med_conf_PELM = merged[merged['_merge'] == 'left_only']
    med_conf_PELM = med_conf_PELM.drop(columns=['_merge'])
    med_conf_PELM['Position_x'] = '-'
    med_conf_PELM.rename(columns={'Position': 'Position_y'}, inplace=True)
    med_conf_PELM = med_conf_PELM.iloc[:, [0, 3, 2, 13, 4, 1]]
    med_conf_PELM['Confidence'] = 'Medium'
    # Convert tuple of string to string
    med_conf_PELM['KinaseID'] = med_conf_PELM['KinaseID'].apply(lambda x: x[0])
    # PhosphoSitePlus
    merged = PhosphoSitePlus_exl.merge(PhosphoSitePlus_i, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'], how='left', indicator=True)
    # Filter rows that are only in PhosphoSitePlus, but not in interactome
    med_conf_PSP = merged[merged['_merge'] == 'left_only']
    med_conf_PSP = med_conf_PSP.drop(columns=['_merge'])
    med_conf_PSP['Position_y'] = '-'
    med_conf_PSP.rename(columns={'Position': 'Position_x'}, inplace=True)
    med_conf_PSP = med_conf_PSP.iloc[:, [0, 3, 2, 1, 4, 9]]
    med_conf_PSP['Confidence'] = 'Medium'
    # Additions
    additions = pd.concat([high_conf, med_conf_PELM, med_conf_PSP])
    PC_prop = pd.concat([PC_prop,
                         pd.DataFrame({'Database': [i],
                                       'Pathway Commons': [(PC_format['INTERACTION_DATA_SOURCE'].str.contains(i)).sum()],
                                       'PhosphoELM_PC': [len(PhosphoELM_i.index)],
                                       'PhosphoSitePlus_PC': [len(PhosphoSitePlus_i.index)],
                                       'High confidence additions': [len(high_conf.index)],
                                       'Medium confidence additions (P.ELM)': [len(med_conf_PELM.index)],
                                       'Medium confidence additions (PSP)': [len(med_conf_PSP.index)]}).set_index("Database")])
    PhosphoELM_i.to_csv(f'..\Results\Pathway_Commons\PhosphoELM_{i}.tsv', sep="\t", index=False)
    PhosphoSitePlus_i.to_csv(f'..\Results\Pathway_Commons\PhosphoSitePlus_{i}.tsv', sep="\t", index=False)
    additions.to_csv(f'..\Results\Additions\{i}_PC.tsv', sep="\t", index=False)
PC_prop.sort_values(by=['Pathway Commons', 'High confidence additions', 'Medium confidence additions (P.ELM)', 'Medium confidence additions (PSP)'], inplace=True, ascending=False)
# PC_prop.to_csv(r'..\Results\Pathway_Commons\PC_prop.tsv', sep="\t")

##Properties of all the interactomes
Interactomes_prop = PC_prop.rename(columns={'Pathway Commons':'Interactomes', 'PhosphoELM_PC':'PhosphoELM', 'PhosphoSitePlus_PC':'PhosphoSitePlus'})
#STRING
merge1 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoELM_STRING, how='inner', on=['SubstrateID','Code','KinaseID','Peptide'])
merge1 = merge1.iloc[:,[0,1,2,3,4,5]]
merge2 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoSitePlus_STRING, how='inner', on=['SubstrateID','Code','KinaseID','Peptide'])
merge2 = merge2.iloc[:, [0, 1, 2, 3, 4, 5]]
overlaps_conf = pd.concat([merge1, merge2])
overlaps_conf.drop_duplicates(inplace=True)
# High confidence additions
merged = PhosphoELMvsPhosphoSitePlus_conf.merge(overlaps_conf, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position_x', 'Position_y'], how='left', indicator=True)
# Filter rows that are only in Phospho.ELM and PhosphoSitePlus, but not in interactome
high_conf = merged[merged['_merge'] == 'left_only']
high_conf = high_conf.drop(columns=['_merge'])
high_conf['Confidence'] = 'High'
high_conf = high_conf.iloc[:, :6].join(high_conf.iloc[:, -1])
# Medium confidence additions
# Phospho.ELM
merged = pd.merge(PhosphoELM_exl, PhosphoELM_STRING, how='inner', on=['SubstrateID', 'Code', 'Peptide','Position'])
merged['KinaseID'] = merged.apply(lambda row: find_overlaps2(row['KinaseID_y'], row['KinaseID_x']),axis=1)
merged.dropna(subset="KinaseID", inplace=True)
merged = merged.iloc[:, :5]
merged.rename(columns={'KinaseID_x': 'KinaseID'}, inplace=True)
med_conf_PELM = PhosphoELM_exl.merge(merged, on=['SubstrateID', 'KinaseID', 'Code', 'Peptide', 'Position'], how='left', indicator=True)
# Filter rows that are only in Phospho.ELM, but not in interactome
med_conf_PELM = med_conf_PELM[med_conf_PELM['_merge'] == 'left_only']
med_conf_PELM = med_conf_PELM.drop(columns=['_merge'])
med_conf_PELM['Position_x'] = '-'
med_conf_PELM.rename(columns={'Position': 'Position_y'}, inplace=True)
med_conf_PELM = med_conf_PELM.iloc[:, [0, 3, 2, 7, 4, 1]]
med_conf_PELM['Confidence'] = 'Medium'
# Convert tuple of string to string
med_conf_PELM['KinaseID'] = med_conf_PELM['KinaseID'].apply(lambda x: x[0])
# PhosphoSitePlus
merged = PhosphoSitePlus_exl.merge(PhosphoSitePlus_STRING, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'], how='left', indicator=True)
# Filter rows that are only in PhosphoSitePlus, but not in interactome
med_conf_PSP = merged[merged['_merge'] == 'left_only']
med_conf_PSP = med_conf_PSP.drop(columns=['_merge'])
med_conf_PSP['Position_y'] = '-'
med_conf_PSP.rename(columns={'Position': 'Position_x'}, inplace=True)
med_conf_PSP = med_conf_PSP.iloc[:, [0, 3, 2, 1, 4, 6]]
med_conf_PSP['Confidence'] = 'Medium'
# Additions
additions = pd.concat([high_conf, med_conf_PELM, med_conf_PSP])
Interactomes_prop = pd.concat([Interactomes_prop,
                               pd.DataFrame({'Database': ['STRING'],
                                             'Interactomes': [len(STRING_format.index)],
                                             'PhosphoELM': [len(PhosphoELM_STRING.index)],
                                             'PhosphoSitePlus': [len(PhosphoSitePlus_STRING.index)],
                                             'High confidence additions': [len(high_conf.index)],
                                             'Medium confidence additions (P.ELM)': [len(med_conf_PELM.index)],
                                             'Medium confidence additions (PSP)': [len(med_conf_PSP.index)]}).set_index("Database")])
additions.to_csv(f'..\Results\Additions\STRING.tsv', sep="\t", index=False)

#BioGRID
merge1 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoELM_BioGRID, how='inner', on=['SubstrateID', 'Code', 'Peptide'])
merge1['KinaseID'] = merge1.apply(lambda row: find_overlaps2(row['KinaseID_x'], row['KinaseID_y']),axis=1)
merge1.dropna(subset="KinaseID", inplace=True)
merge1 = merge1.iloc[:, [0, 16, 2, 3, 4, 5]]
merge2 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoSitePlus_BioGRID, how='inner', on=['SubstrateID','Code','KinaseID','Peptide'])
merge2 = merge2.iloc[:, [0, 1, 2, 3, 4, 5]]
overlaps_conf = pd.concat([merge1, merge2])
overlaps_conf.drop_duplicates(inplace=True)
# High confidence additions
merged = PhosphoELMvsPhosphoSitePlus_conf.merge(overlaps_conf, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position_x', 'Position_y'], how='left', indicator=True)
# Filter rows that are only in Phospho.ELM and PhosphoSitePlus, but not in interactome
high_conf = merged[merged['_merge'] == 'left_only']
high_conf = high_conf.drop(columns=['_merge'])
high_conf['Confidence'] = 'High'
high_conf = high_conf.iloc[:, :6].join(high_conf.iloc[:, -1])
# Medium confidence additions
# Phospho.ELM
merged = pd.merge(PhosphoELM_exl, PhosphoELM_BioGRID, how='inner', on=['SubstrateID', 'Code', 'Peptide','Position'])
merged['KinaseID'] = merged.apply(lambda row: find_overlaps(row['KinaseID_x'], row['KinaseID_y']), axis=1)
merged = merged[merged['KinaseID'].apply(lambda x: len(x) > 0)]
merged = merged.iloc[:, :5]
merged.rename(columns={'KinaseID_x': 'KinaseID'}, inplace=True)
med_conf_PELM = PhosphoELM_exl.merge(merged, on=['SubstrateID', 'KinaseID', 'Code', 'Peptide', 'Position'], how='left', indicator=True)
# Filter rows that are only in Phospho.ELM, but not in interactome
med_conf_PELM = med_conf_PELM[med_conf_PELM['_merge'] == 'left_only']
med_conf_PELM = med_conf_PELM.drop(columns=['_merge'])
med_conf_PELM['Position_x'] = '-'
med_conf_PELM.rename(columns={'Position': 'Position_y'}, inplace=True)
med_conf_PELM = med_conf_PELM.iloc[:, [0, 3, 2, 7, 4, 1]]
med_conf_PELM['Confidence'] = 'Medium'
# Convert tuple of string to string
med_conf_PELM['KinaseID'] = med_conf_PELM['KinaseID'].apply(lambda x: x[0])
# PhosphoSitePlus
merged = PhosphoSitePlus_exl.merge(PhosphoSitePlus_BioGRID, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'], how='left', indicator=True)
# Filter rows that are only in PhosphoSitePlus, but not in interactome
med_conf_PSP = merged[merged['_merge'] == 'left_only']
med_conf_PSP = med_conf_PSP.drop(columns=['_merge'])
med_conf_PSP['Position_y'] = '-'
med_conf_PSP.rename(columns={'Position': 'Position_x'}, inplace=True)
med_conf_PSP = med_conf_PSP.iloc[:, [0, 3, 2, 1, 4, 7]]
med_conf_PSP['Confidence'] = 'Medium'
# Additions
additions = pd.concat([high_conf, med_conf_PELM, med_conf_PSP])
Interactomes_prop = pd.concat([Interactomes_prop,
                               pd.DataFrame({'Database': ['BioGRID'],
                                             'Interactomes': [len(BioGRID_format.index)],
                                             'PhosphoELM': [len(PhosphoELM_BioGRID.index)],
                                             'PhosphoSitePlus': [len(PhosphoSitePlus_BioGRID.index)],
                                             'High confidence additions': [len(high_conf.index)],
                                             'Medium confidence additions (P.ELM)': [len(med_conf_PELM.index)],
                                             'Medium confidence additions (PSP)': [len(med_conf_PSP.index)]}).set_index("Database")])
additions.to_csv(f'..\Results\Additions\BioGRID.tsv', sep="\t", index=False)
Interactomes_prop["Total additions"] = Interactomes_prop['High confidence additions'] + Interactomes_prop['Medium confidence additions (P.ELM)'] + Interactomes_prop['Medium confidence additions (PSP)']
Interactomes_prop.sort_values(by=['Interactomes', 'Total additions'], inplace=True, ascending=False)
# Interactomes_prop.to_csv(r'..\Results\Interactomes_prop.tsv', sep="\t")

##Bar graph for all the interactomes
matplotlib.use('Qt5Agg') #To prevent matplotlib from not responding
max_val = Interactomes_prop['Total additions'].max()
y_values = range(0, max_val+2000, 1000)
#green:STRING, red:BioGRID and blue:Pathway Commons
colors = ['green'] + ['red'] + ['blue']*14
fig, ax = plt.subplots(figsize=(8, 6))
patterns = ['', '.', '//']
hatch_legend_labels = ['High', 'Medium (Phospho.ELM)', 'Medium (PhosphoSitePlus)']
for i, row in enumerate(Interactomes_prop.iterrows()):
    ax.bar(i, row[1]['Medium confidence additions (PSP)'], color=colors[i], hatch='//', alpha =0.6)
    ax.bar(i, row[1]['Medium confidence additions (P.ELM)'], color=colors[i], hatch='.', alpha = 0.6, bottom=row[1]['Medium confidence additions (PSP)'])
    ax.bar(i, row[1]['High confidence additions'], color=colors[i], bottom=row[1]['Medium confidence additions (PSP)']+row[1]['Medium confidence additions (P.ELM)'])
ax.set_xticks(range(len(Interactomes_prop)))
ax.set_xticklabels(Interactomes_prop.index, rotation=45)
ax.set_yticks(y_values)
ax.set_yticklabels(y_values)
ax.tick_params(axis='x', labelsize=8)
plt.ylabel('Number of Additions')
green_patch = plt.Rectangle((0,0), 1, 1, color='green')
red_patch = plt.Rectangle((0,0), 1, 1, color='red')
blue_patch = plt.Rectangle((0,0), 1, 1, color='blue')
legend1 = ax.legend([green_patch,red_patch,blue_patch], ['STRING','BioGRID','Pathway Commons'],bbox_to_anchor=(0.8, 0.98))
legend_handles = [Patch(facecolor='white', edgecolor='grey', hatch=patterns[i]) for i in range(len(patterns))]
legend2 = ax.legend(legend_handles, hatch_legend_labels,title='Confidence')
# Add both legends to the plot
ax.add_artist(legend1)
plt.savefig('..\Results\Additions\Interactomes_additions.png')

##Venn diagrams for the visualisations of the overlaps
#Labels are the number of edges (predictions)
#Phospho.ELM vs STRING
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(STRING_format.index)*0.1, len(PhosphoELM_conf.index), len(PhosphoELM_STRING.index)),
             set_colors = ('green', 'brown'),
             set_labels = ('STRING', 'Phospho.ELM'))
ax.set_title('Phospho.ELM vs STRING')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['STRING','Edges']-len(PhosphoELM_STRING.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(PhosphoELM_STRING.index)}")
venn.get_label_by_id('11').set_text(f"{len(PhosphoELM_STRING.index)}")
venn.get_label_by_id('01').set_position((0.61,0))
plt.legend([len(STRING_format.index), len(PhosphoELM_conf.index)],title="Total interactions")
plt.savefig('..\Results\Figures\Phospho.ELMvsSTRING.png')

#Phospho.ELM vs BioGRID
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(BioGRID_format.index), len(PhosphoELM_conf.index), len(PhosphoELM_BioGRID.index)),
             set_colors = ('red', 'brown'),
             set_labels = ('BioGRID', 'Phospho.ELM'))
ax.set_title('Phospho.ELM vs BioGRID')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['BioGRID','Edges']-len(PhosphoELM_BioGRID.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['Phospho.ELM','Edges']-len(PhosphoELM_BioGRID.index)}")
venn.get_label_by_id('11').set_text(f"{len(PhosphoELM_BioGRID.index)}")
plt.legend([len(BioGRID_format.index), len(PhosphoELM_conf.index)],title="Total interactions")
plt.savefig('..\Results\Figures\Phospho.ELMvsBioGRID.png')

#Phospho.ELM vs Pathway Commons
PC_total = PC_format[~PC_format.apply(lambda x: tuple(sorted(x[['protein1', 'protein2']])), axis=1).duplicated(keep='first')]

PhosphoELM_PC = PhosphoELM_PC.drop(['INTERACTION_DATA_SOURCE','INTERACTION_TYPE','PARTICIPANT_A','PARTICIPANT_B'], axis=1)
PhosphoELM_PC.drop_duplicates(inplace=True)
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PC_total.index), len(PhosphoELM_conf.index), len(PhosphoELM_PC.index)),
             set_colors = ('grey', 'brown'),
             set_labels = ('Pathway Commons', 'Phospho.ELM'))
ax.set_title('Phospho.ELM vs Pathway Commons')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Pathway Commons','Edges']-len(PhosphoELM_PC.index)}")
venn.get_label_by_id('01').set_text(f"{len(PhosphoELM_conf.index)-len(PhosphoELM_PC.index)}")
venn.get_label_by_id('11').set_text(f"{len(PhosphoELM_PC.index)}")
venn.get_label_by_id('01').set_position((0.63,0))
plt.legend([len(PC_total.index), len(PhosphoELM_conf.index)],title="Total interactions")
plt.savefig(r'..\Results\Figures\Phospho.ELMvsPC.png')

#PhosphoSitePlus vs STRING
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(STRING_format.index)*0.1, len(PhosphoSitePlus_format.index), len(PhosphoSitePlus_STRING.index)),
             set_colors = ('green', 'orange'),
             set_labels = ('STRING', 'PhosphoSitePlus'))
ax.set_title('PhosphoSitePlus vs STRING')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['STRING','Edges']-len(PhosphoSitePlus_STRING.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(PhosphoSitePlus_STRING.index)}")
venn.get_label_by_id('11').set_text(f"{len(PhosphoSitePlus_STRING.index)}")
plt.legend([len(STRING_format.index), len(PhosphoSitePlus_format.index)],title="Total interactions")
plt.savefig('..\Results\Figures\PhosphoSitePlusvsSTRING.png')

#PhosphoSitePlus vs BioGRID
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(BioGRID_format.index), len(PhosphoSitePlus_format.index), len(PhosphoSitePlus_BioGRID.index)),
             set_colors = ('red', 'orange'),
             set_labels = ('BioGRID', 'PhosphoSitePlus'))
ax.set_title('PhosphoSitePlus vs BioGRID')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['BioGRID','Edges']-len(PhosphoSitePlus_BioGRID.index)}")
venn.get_label_by_id('01').set_text(f"{net_prop.loc['PhosphoSitePlus','Edges']-len(PhosphoSitePlus_BioGRID.index)}")
venn.get_label_by_id('11').set_text(f"{len(PhosphoSitePlus_BioGRID.index)}")
plt.legend([len(BioGRID_format.index), len(PhosphoSitePlus_format.index)],title="Total interactions")
plt.savefig('..\Results\Figures\PhosphoSitePlusvsBioGRID.png')

#PhosphoSitePlus vs Pathway Commons
PhosphoSitePlus_PC = PhosphoSitePlus_PC.drop(['INTERACTION_DATA_SOURCE','INTERACTION_TYPE','PARTICIPANT_A','PARTICIPANT_B'], axis=1)
PhosphoSitePlus_PC.drop_duplicates(inplace=True)
fig, ax = plt.subplots(figsize=(8, 6))
venn = venn2(subsets = (len(PC_total.index), len(PhosphoSitePlus_format.index), len(PhosphoSitePlus_PC.index)),
             set_colors = ('grey', 'orange'),
             set_labels = ('Pathway Commons', 'PhosphoSitePlus'))
ax.set_title('PhosphoSitePlus vs Pathway Commons')
venn.get_patch_by_id('11').set_color('yellow')
venn.get_label_by_id('10').set_text(f"{net_prop.loc['Pathway Commons','Edges']-len(PhosphoSitePlus_PC.index)}")
venn.get_label_by_id('01').set_text(f"{len(PhosphoSitePlus_format.index)-len(PhosphoSitePlus_PC.index)}")
venn.get_label_by_id('11').set_text(f"{len(PhosphoSitePlus_PC.index)}")
plt.legend([len(PC_total.index), len(PhosphoSitePlus_format.index)],title="Total interactions")
plt.savefig(r'..\Results\Figures\PhosphoSitePlusvsPC.png')

##Classify interaction types for each of the individual databases of Pathway Commons
interaction_types=pd.DataFrame()
itypes = PC_format['INTERACTION_TYPE'].unique()
for i in PC_databases:
    PC_i = PC_format[PC_format['INTERACTION_DATA_SOURCE'].str.contains(i)]
    interaction_types = pd.concat([interaction_types,
                         pd.DataFrame({'Database': [i],
                                       'Pathway Commons': [len(PC_i.index)],
                                       itypes[0]: [(PC_i['INTERACTION_TYPE'].str.contains(itypes[0])).sum()],
                                       itypes[1]: [(PC_i['INTERACTION_TYPE'].str.contains(itypes[1])).sum()],
                                       itypes[2]: [(PC_i['INTERACTION_TYPE'].str.contains(itypes[2])).sum()],
                                       itypes[3]: [(PC_i['INTERACTION_TYPE'].str.contains(itypes[3])).sum()],
                                       itypes[4]: [(PC_i['INTERACTION_TYPE'].str.contains(itypes[4])).sum()]}).set_index("Database")])
interaction_types.sort_values(by=['Pathway Commons'], inplace=True, ascending=False)
# interaction_types.to_csv(r'..\Results\Pathway_Commons\interaction_types.tsv', sep="\t")

##Bar graph for all the interaction types of Pathway Commons
max_val = interaction_types['Pathway Commons'].max()
y_values = range(0, max_val+1000, 5000)
ax = interaction_types.iloc[:, -5:].plot(kind='bar', stacked=True, figsize=(8, 6))
ax.set_xticks(range(len(interaction_types)))
ax.set_xticklabels(interaction_types.index, rotation=45)
ax.set_yticks(y_values)
ax.set_yticklabels(y_values)
ax.tick_params(axis='x', labelsize=8)
plt.ylabel('Number of interactions')
plt.legend()
plt.savefig('..\Results\Pathway_Commons\interactions_types.png')

# Calculate the percentage of each column with respect to the total number of interactions
df_percentage = interaction_types.div(interaction_types['Pathway Commons'], axis=0) * 100
y_values = range(0, 110, 10)
ax = df_percentage.iloc[:, -5:].plot(kind='bar', stacked=True, figsize=(10, 9))
ax.set_xticks(range(len(interaction_types)))
ax.set_xticklabels(interaction_types.index, rotation=45)
ax.set_yticks(y_values)
ax.set_yticklabels(y_values)
ax.tick_params(axis='x', labelsize=8)
plt.ylabel('Percentage')
plt.legend(bbox_to_anchor=(0.8, 0.98))
plt.savefig('..\Results\Pathway_Commons\interactions_types_percentage.png')

## Represent Interactomes properties as percentages (STRING, BioGRID and Pathway Commons)
# Get gene names of proteins in STRING
# Map UniProt IDs to gene names
# region uniprot_to_gname
# uniprot = list(set(STRING_format['protein1'].tolist()+STRING_format['protein2'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# # uniprot_to_gname = {}
# for uni in uniprot[16000:18000]:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{uni}&fields=gene_primary'
#     response = requests.get(query_url)
#     gene = response.text.split('\n')[1:-1]
#     uniprot_to_gname[uni] = gene
# with open(r'.\Variables\uniprot_to_gname.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_gname, f)
# endregion
STRING = pd.DataFrame()
STRING['protein1'] = STRING_format.apply(lambda x: uniprot_to_gname[x['protein1']][0] if uniprot_to_gname[x['protein1']]!=[] else pd.NA, axis=1)
STRING['protein2'] = STRING_format.apply(lambda x: uniprot_to_gname[x['protein2']][0] if uniprot_to_gname[x['protein2']]!=[] else pd.NA, axis=1)
STRING.dropna(inplace=True)
#103 interactions are dropped
BioGRID = BioGRID_format.iloc[:,:2]
BioGRID = BioGRID.copy()
BioGRID.columns=['protein1', 'protein2']
PC = PC_format.iloc[:, :2]
PC = PC.copy()
PC.columns=['protein1', 'protein2']
# Total PPIs in the interactomes
total_ppi = pd.concat([STRING, BioGRID, PC])
# total_ppi = total_ppi[~total_ppi.apply(lambda x: tuple(sorted(x[['protein1', 'protein2']])), axis=1).duplicated(keep='first')]
total_ppi.to_csv(r'..\Results\total_ppi.tsv', sep="\t", index=False)

#Exclusive interactions
#STRING
BioGRID_PC = pd.concat([BioGRID, PC])
BioGRID_PC.drop_duplicates(inplace=True)
merge1 = pd.merge(STRING, BioGRID_PC, how='left', indicator=True)
STRING_reverse = STRING.copy()
STRING_reverse.columns=['protein2', 'protein1']
merge2 = pd.merge(STRING_reverse, BioGRID_PC, how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only']
merge1 = merge1.drop(columns=['_merge'])
merge2 = merge2[merge2['_merge'] == 'left_only']
merge2 = merge2.drop(columns=['_merge'])
merge2.columns=['protein1', 'protein2']
STRING_excl = pd.concat([merge1,merge2])
STRING_excl.drop_duplicates(inplace=True)
STRING_excl.to_csv(r'..\Results\STRING_exclusive_ppi.tsv', sep="\t", index=False)

#BioGRID
PC_exclude_BioGRID = PC_format[PC_format['INTERACTION_DATA_SOURCE']!='BioGRID']
PC_exclude_BioGRID = PC_exclude_BioGRID.iloc[:,:2]
PC_exclude_BioGRID.columns=['protein1', 'protein2']
STRING_PC = pd.concat([STRING, PC_exclude_BioGRID])
STRING_PC.drop_duplicates(inplace=True)
merge1 = pd.merge(BioGRID, STRING_PC, how='left', indicator=True)
BioGRID_reverse = BioGRID.copy()
BioGRID_reverse.columns=['protein2', 'protein1']
merge2 = pd.merge(BioGRID_reverse, STRING_PC, how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only']
merge1 = merge1.drop(columns=['_merge'])
merge2 = merge2[merge2['_merge'] == 'left_only']
merge2 = merge2.drop(columns=['_merge'])
merge2.columns=['protein1', 'protein2']
BioGRID_excl = pd.concat([merge1,merge2])
BioGRID_excl.drop_duplicates(inplace=True)
BioGRID_excl.to_csv(r'..\Results\BioGRID_exclusive_ppi.tsv', sep="\t", index=False)

#PC
STRING_BioGRID = pd.concat([STRING, BioGRID])
STRING_BioGRID.drop_duplicates(inplace=True)
merge1 = pd.merge(PC, STRING_BioGRID, how='left', indicator=True)
PC_reverse = PC.copy()
PC_reverse.columns=['protein2', 'protein1']
merge2 = pd.merge(PC_reverse, STRING_BioGRID, how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only']
merge1 = merge1.drop(columns=['_merge'])
merge2 = merge2[merge2['_merge'] == 'left_only']
merge2 = merge2.drop(columns=['_merge'])
merge2.columns=['protein1', 'protein2']
PC_excl = pd.concat([merge1,merge2])
PC_excl.drop_duplicates(inplace=True)
PC_excl.to_csv(r'..\Results\PC_exclusive_ppi.tsv', sep="\t", index=False)

# Get statistics for PC versus Phospho.ELM and PhosphoSitePlus
# Retain only unique interactions
PhosphoELM_i = PhosphoELM_PC.iloc[:,:9].drop_duplicates(subset=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'])
PhosphoSitePlus_i = PhosphoSitePlus_PC.iloc[:,:7].drop_duplicates(subset=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'])
merge1 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoELM_i, how='inner', on=['SubstrateID', 'Code', 'Peptide'])
if not merge1.empty:
    merge1['KinaseID'] = merge1.apply(lambda row: find_overlaps2(row['KinaseID_x'], row['KinaseID_y']), axis=1)
    merge1.dropna(subset="KinaseID", inplace=True)
    merge1 = merge1.iloc[:, [0, 16, 2, 3, 4, 5]]
else:
    merge1 = pd.DataFrame()
merge2 = pd.merge(PhosphoELMvsPhosphoSitePlus_conf, PhosphoSitePlus_i, how='inner',
                  on=['SubstrateID', 'Code', 'KinaseID', 'Peptide'])
merge2 = merge2.iloc[:, [0, 1, 2, 3, 4, 5]]
overlaps_conf = pd.concat([merge1, merge2])
overlaps_conf.drop_duplicates(inplace=True)
# High confidence additions
merged = PhosphoELMvsPhosphoSitePlus_conf.merge(overlaps_conf,
                                                on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position_x',
                                                    'Position_y'], how='left', indicator=True)
# Filter rows that are only in Phospho.ELM and PhosphoSitePlus, but not in interactome
high_conf = merged[merged['_merge'] == 'left_only']
high_conf = high_conf.drop(columns=['_merge'])
high_conf['Confidence'] = 'High'
high_conf = high_conf.iloc[:, :6].join(high_conf.iloc[:, -1])
# Medium confidence additions
# Phospho.ELM
merged = PhosphoELM_exl.merge(PhosphoELM_i, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'], how='left',
                              indicator=True)
# Filter rows that are only in Phospho.ELM, but not in interactome
med_conf_PELM = merged[merged['_merge'] == 'left_only']
med_conf_PELM = med_conf_PELM.drop(columns=['_merge'])
med_conf_PELM['Position_x'] = '-'
med_conf_PELM.rename(columns={'Position': 'Position_y'}, inplace=True)
med_conf_PELM = med_conf_PELM.iloc[:, [0, 3, 2, 11, 4, 1]]
med_conf_PELM['Confidence'] = 'Medium'
# Convert tuple of string to string
med_conf_PELM['KinaseID'] = med_conf_PELM['KinaseID'].apply(lambda x: x[0])
# PhosphoSitePlus
merged = PhosphoSitePlus_exl.merge(PhosphoSitePlus_i, on=['SubstrateID', 'Code', 'KinaseID', 'Peptide', 'Position'],
                                   how='left', indicator=True)
# Filter rows that are only in PhosphoSitePlus, but not in interactome
med_conf_PSP = merged[merged['_merge'] == 'left_only']
med_conf_PSP = med_conf_PSP.drop(columns=['_merge'])
med_conf_PSP['Position_y'] = '-'
med_conf_PSP.rename(columns={'Position': 'Position_x'}, inplace=True)
med_conf_PSP = med_conf_PSP.iloc[:, [0, 3, 2, 1, 4, 7]]
med_conf_PSP['Confidence'] = 'Medium'
# Additions
additions = pd.concat([high_conf, med_conf_PELM, med_conf_PSP])
additions.to_csv(f'..\Results\Additions\Pathway_Commons.tsv', sep="\t", index=False)
#Autophosphorylation interactions in additions
STRING_additions = pd.read_csv(f'..\Results\Additions\STRING.tsv', sep="\t")
BioGRID_additions = pd.read_csv(f'..\Results\Additions\BioGRID.tsv', sep="\t")
PC_additions = pd.read_csv(f'..\Results\Additions\Pathway_Commons.tsv', sep="\t")

Interactomes_prop_per = pd.DataFrame()
Interactomes_prop_per = pd.concat([Interactomes_prop_per,
                         pd.DataFrame({'Database': ['STRING','BioGRID','Pathway Commons'],
                                       '# K-S interactions': [len(STRING.index), len(BioGRID.index), len(PC.index)],
                                       '% coverage of K-S interactions': [round(len(STRING.index)/len(total_ppi.index)*100,1), round(len(BioGRID.index)/len(total_ppi.index)*100,1), round(len(PC.index)/len(total_ppi.index)*100,1)],
                                       '% of exclusive K-S interactions': [round(len(STRING_excl.index)/len(total_ppi.index)*100,1), round(len(BioGRID_excl.index)/len(total_ppi.index)*100,1), round(len(PC_excl.index)/len(total_ppi.index)*100,1)],
                                       'PhosphoELM_PC': [Interactomes_prop.iloc[0,1], Interactomes_prop.iloc[1,1], len(PhosphoELM_i.index)],
                                       'PhosphoSitePlus_PC': [Interactomes_prop.iloc[0,2], Interactomes_prop.iloc[1,2], len(PhosphoSitePlus_i.index)],
                                       'High confidence additions': [Interactomes_prop.iloc[0,3], Interactomes_prop.iloc[1,3], len(high_conf.index)],
                                       'Medium confidence additions (P.ELM)': [Interactomes_prop.iloc[0,4], Interactomes_prop.iloc[1,4], len(med_conf_PELM.index)],
                                       'Medium confidence additions (PSP)': [Interactomes_prop.iloc[0,5], Interactomes_prop.iloc[1,5], len(med_conf_PSP.index)],
                                       'Total additions': [Interactomes_prop.iloc[0,6], Interactomes_prop.iloc[1,6], len(additions.index)],
                                        'Autophosphorylation': [len(STRING_additions[STRING_additions['SubstrateID']==STRING_additions['KinaseID']]), len(BioGRID_additions[BioGRID_additions['SubstrateID']==BioGRID_additions['KinaseID']]), len(PC_additions[PC_additions['SubstrateID']==PC_additions['KinaseID']])]}).set_index("Database")])
Interactomes_prop_per.sort_values(by=['# K-S interactions', 'Total additions'], inplace=True, ascending=False)
# Interactomes_prop_per.to_csv(r'..\Results\Interactomes_prop_per.tsv', sep="\t")

##Bar graph for all the interactomes
y_values = range(0, 101, 10)
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticks(range(len(Interactomes_prop_per)))
ax.set_yticks(y_values)
ax.set_yticklabels(y_values)
bar_plot = Interactomes_prop_per.plot(y=["% coverage of K-S interactions", "% of exclusive K-S interactions"], kind="bar", ax=ax)
for p in bar_plot.patches:
    ax.annotate(f'{p.get_height()}', (p.get_x() + p.get_width() / 2, p.get_height()),ha='center', va='center', xytext=(0, 10), textcoords='offset points')
ax.legend(["Total K-S interactions", "Exclusive K-S interactions"])
ax.set_xticklabels(Interactomes_prop_per.index, rotation=360)
ax.set_xlabel('')
ax.set_ylim(0, 100)
plt.ylabel('Percentage of PPIs')
plt.savefig('..\Results\Figures\Interactomes_per_ppi.png')

##Bar graph for confidence in the additional interactions in all the interactomes
y_values = range(0, 101, 10)
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticks(range(len(Interactomes_prop_per)))
ax.set_yticks(y_values)
ax.set_yticklabels(y_values)
df = Interactomes_prop_per.copy()
df["High confidence additions (%)"] = round((Interactomes_prop_per["High confidence additions"] / Interactomes_prop_per["Total additions"]) * 100,1)
df["Medium confidence additions (P.ELM) (%)"] = round((Interactomes_prop_per["Medium confidence additions (P.ELM)"] / Interactomes_prop_per["Total additions"]) * 100,1)
df["Medium confidence additions (PSP) (%)"] = round((Interactomes_prop_per["Medium confidence additions (PSP)"] / Interactomes_prop_per["Total additions"]) * 100,1)
bar_plot = df.plot(y=["High confidence additions (%)", "Medium confidence additions (P.ELM) (%)", "Medium confidence additions (PSP) (%)"], kind="bar", ax=ax)
ax.legend(["High", "Medium (Phospho.ELM)", "Medium (PhosphoSitePlus)"], title='Confidence', bbox_to_anchor=(0.73, 0.94))
for p in bar_plot.patches:
    ax.annotate(f'{p.get_height()}', (p.get_x() + p.get_width() / 2, p.get_height()),ha='center', va='center', xytext=(0, 5), textcoords='offset points')
ax.set_xticklabels(Interactomes_prop_per.index, rotation=360)
ax.set_xlabel('')
ax.set_ylim(0, 100)
plt.ylabel('Percentage of additions')
plt.savefig('..\Results\Additions\Interactomes_additions_per.png')

## Interactions for network analysis
## Include the additional interactions of BioGRID (2023) to Pathway Commons
overlaps2 = overlaps.iloc[:,:4]
BioGRID_additional = pd.merge(BioGRID_format, overlaps2, how='left', indicator=True)
BioGRID_additional = BioGRID_additional[BioGRID_additional['_merge'] == 'left_only'].drop(columns=['_merge'])
BioGRID_additional = BioGRID_additional.iloc[:,:2]
BioGRID_additional.columns=['protein1', 'protein2']
PC_new = pd.concat([PC, BioGRID_additional])
PC_new.drop_duplicates(inplace=True)
# Exclusive interactions
merge1 = pd.merge(PC_new, STRING, how='left', indicator=True)
PC_new_reverse = PC_new.copy()
PC_new_reverse.columns=['protein2', 'protein1']
merge2 = pd.merge(PC_new_reverse, STRING, how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2.columns=['protein1', 'protein2']
PC_new_BioGRID_excl = pd.concat([merge1,merge2])
PC_new_BioGRID_excl.drop_duplicates(inplace=True)
PC_new_BioGRID_excl.to_csv(r'..\Results\Network_Analysis\PC_BioGRID_exclusive_ppi.tsv', sep="\t", index=False)
# Get statistics for PC_new (PC + BioGRID) versus Phospho.ELM and PhosphoSitePlus
PC_additions = pd.read_csv(f'..\Results\Additions\Pathway_Commons.tsv', sep="\t")
BioGRID_additions = pd.read_csv(f'..\Results\Additions\BioGRID.tsv', sep="\t")
additions = pd.concat([PC_additions, BioGRID_additions])
additions.drop_duplicates(inplace=True)
additions.sort_values(by=['Confidence', 'Position_x', 'Position_y'], inplace=True)
additions.to_csv(r'..\Results\Network_Analysis\PC_BioGRID_additions.tsv', sep="\t", index=False)

## Include the BioGRID (2023) to Reactome alone
Reactome = PC_format[PC_format['INTERACTION_DATA_SOURCE'].str.contains("Reactome")]
BioGRID_Reactome = pd.concat([BioGRID, Reactome.iloc[:,[2,3]]])
# Exclusive interactions
Reactome = Reactome[Reactome['INTERACTION_DATA_SOURCE']=='Reactome']
Reactome = Reactome.iloc[:,:2]
Reactome.columns=['protein1', 'protein2']
merge1 = pd.merge(Reactome, STRING, how='left', indicator=True)
Reactome_reverse = Reactome.copy()
Reactome_reverse.columns=['protein2', 'protein1']
merge2 = pd.merge(Reactome_reverse, STRING, how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2.columns=['protein1', 'protein2']
Reactome_BioGRID_excl = pd.concat([merge1,merge2,BioGRID_excl])
Reactome_BioGRID_excl.drop_duplicates(inplace=True)
Reactome_BioGRID_excl.to_csv(r'..\Results\Network_Analysis\Reactome_exclusive_ppi.tsv', sep="\t", index=False)
# Get statistics for (Reactome + BioGRID) versus Phospho.ELM and PhosphoSitePlus
Reactome_additions = pd.read_csv(f'..\Results\Additions\Reactome_PC.tsv', sep="\t")
additions1 = pd.concat([Reactome_additions, BioGRID_additions])
additions1.drop_duplicates(inplace=True)
additions1.sort_values(by=['Confidence', 'Position_x', 'Position_y'], inplace=True)
additions1.to_csv(r'..\Results\Network_Analysis\Reactome_BioGRID_additions.tsv', sep="\t", index=False)

NA_prop_per = pd.DataFrame()
NA_prop_per = pd.concat([NA_prop_per,
                         pd.DataFrame({'Interactomes': ['Pathway Commons + BioGRID', 'BioGRID + Reactome'],
                                       '# K-S interactions': [len(PC_new.index), len(BioGRID_Reactome.index)],
                                       '% coverage of K-S interactions': [round(len(PC_new.index)/len(total_ppi.index)*100,1), round(len(BioGRID_Reactome.index)/len(total_ppi.index)*100,1)],
                                       '% of exclusive K-S interactions': [round(len(PC_new_BioGRID_excl.index)/len(total_ppi.index)*100,1), round(len(Reactome_BioGRID_excl.index)/len(total_ppi.index)*100,1)],
                                       'PhosphoELM_PC': [len(additions[additions['Position_y']!='-']), len(additions1[additions1['Position_y']!='-'])],
                                       'PhosphoSitePlus_PC': [len(additions[additions['Position_x']!='-']), len(additions1[additions1['Position_x']!='-'])],
                                       'High confidence additions': [len(additions[additions['Confidence']=='High']), len(additions1[additions1['Confidence']=='High'])],
                                       'Medium confidence additions (P.ELM)': [len(additions[additions['Position_x']=='-']), len(additions1[additions1['Position_x']=='-'])],
                                       'Medium confidence additions (PSP)': [len(additions[additions['Position_y']=='-']), len(additions1[additions1['Position_y']=='-'])],
                                       'Total additions': [len(additions.index), len(additions1.index)],
                                       'Autophosphorylation': [len(additions[additions['SubstrateID']==additions['KinaseID']]), len(additions1[additions1['SubstrateID']==additions1['KinaseID']])]}).set_index("Interactomes")])
# NA_prop_per.to_csv(r'..\Results\Network_Analysis\NA_prop_per.tsv', sep="\t")