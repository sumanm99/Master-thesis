import pandas as pd
import networkx as nx
import pickle
import requests

#Load saved variables
with open(r'.\Variables\gname_to_pfam.pkl', 'rb') as f:
    gname_to_pfam = pickle.load(f)
with open(r'.\Variables\uniprot_to_gname.pkl', 'rb') as f:
    uniprot_to_gname = pickle.load(f)
PC_new_BioGRID_excl = pd.read_csv(r'..\Results\Network_Analysis\PC_BioGRID_exclusive_ppi.tsv', sep="\t")
Reactome_BioGRID_excl = pd.read_csv(r'..\Results\Network_Analysis\Reactome_exclusive_ppi.tsv', sep="\t")
PC_BioGRID_additions = pd.read_csv(r'..\Results\Network_Analysis\PC_BioGRID_additions.tsv', sep="\t")
Reactome_BioGRID_additions = pd.read_csv(r'..\Results\Network_Analysis\Reactome_BioGRID_additions.tsv', sep="\t")

## Breast cancer paper
# Load the breast cancer PPI file
graph = nx.read_graphml("..\Data\Breast cancer protein-protein interaction network.graphml")
edges_data = graph.edges(data=True)
breast_cancer = pd.DataFrame(edges_data, columns=['source', 'target', 'attributes'])
# The key 'name' of the dictionary in the 'attributes' column contains information of the interacting nodes
breast_cancer['name'] = breast_cancer['attributes'].apply(lambda x: x.get('name'))
breast_cancer[['source', 'target']] = breast_cancer['name'].str.split(' \(interacts with\) ', expand=True)
breast_cancer = breast_cancer.iloc[:, :2]
breast_cancer.drop_duplicates(inplace=True)
# print(pd.concat([breast_cancer['source'],breast_cancer['target']]).nunique()) # 523 proteins
#Map gene names to protein families to filter out those which have "kinase"
# region gname_to_pfam
# gene_names= list(set(breast_cancer['source'].tolist()+breast_cancer['target'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# gene_names = [x for x in gene_names if x not in gname_to_pfam.keys()]
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_name = response.text.split('\n')[1]
#     if "kinase" in protein_name.lower():
#         gname_to_pfam[gene_name] = protein_name
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
#Retain phosphorylation interactions
kinases = gname_to_pfam.keys()
breast_cancer_source = breast_cancer[breast_cancer['source'].isin(kinases)]
breast_cancer_target = breast_cancer[breast_cancer['target'].isin(kinases)]
breast_cancer = pd.concat([breast_cancer_source, breast_cancer_target])
breast_cancer.drop_duplicates(inplace=True)

# Check for interactions that are not in the chosen interactome networks
# PC + BioGRID
merge1 = pd.merge(breast_cancer, PC_new_BioGRID_excl, left_on=['source', 'target'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(breast_cancer, PC_new_BioGRID_excl, left_on=['source', 'target'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
breast_cancer_excl = pd.concat([merge1, merge2])
breast_cancer_excl.drop_duplicates(inplace=True)
breast_cancer_excl = breast_cancer_excl.iloc[:,:2]
# Additions that are contributed by the breast cancer data exclusively
# Get gene names of proteins in the additions file
# Map UniProt IDs to gene names
# region uniprot_to_gname
# uniprot = list(set(PC_BioGRID_additions['SubstrateID'].tolist()+PC_BioGRID_additions['KinaseID'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# uniprot = [x for x in uniprot if x not in uniprot_to_gname.keys()]
# for uni in uniprot:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{uni}&fields=gene_primary'
#     response = requests.get(query_url)
#     gene = response.text.split('\n')[1:-1]
#     uniprot_to_gname[uni] = gene
# with open(r'.\Variables\uniprot_to_gname.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_gname, f)
# endregion
PC_BioGRID_additions['protein1'] = PC_BioGRID_additions.apply(lambda x: uniprot_to_gname[x['SubstrateID']][0] if uniprot_to_gname[x['SubstrateID']]!=[] else pd.NA, axis=1)
PC_BioGRID_additions['protein2'] = PC_BioGRID_additions.apply(lambda x: uniprot_to_gname[x['KinaseID']][0] if uniprot_to_gname[x['KinaseID']]!=[] else pd.NA, axis=1)
PC_BioGRID_additions.dropna(inplace=True)
PC_BioGRID_additions = PC_BioGRID_additions.iloc[:,[7,8]]
merge1 = pd.merge(breast_cancer_excl, PC_BioGRID_additions, left_on=['source', 'target'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(breast_cancer_excl, PC_BioGRID_additions, left_on=['source', 'target'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
breast_cancer_additions1 = pd.concat([merge1, merge2])
breast_cancer_additions1.drop_duplicates(inplace=True)
breast_cancer_additions1 = breast_cancer_additions1.iloc[:, :2]
breast_cancer_additions1.to_csv(r'..\Results\Network_Analysis\PC_BioGRID_additions(breast_cancer).tsv', sep="\t", index=False)

# BioGRID + Reactome
merge1 = pd.merge(breast_cancer, Reactome_BioGRID_excl, left_on=['source', 'target'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(breast_cancer, Reactome_BioGRID_excl, left_on=['source', 'target'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
breast_cancer_excl1 = pd.concat([merge1, merge2])
breast_cancer_excl1.drop_duplicates(inplace=True)
breast_cancer_excl1 = breast_cancer_excl1.iloc[:,:2]
# Additions that are contributed by the breast cancer data exclusively
# Get gene names of proteins in the additions file
# Map UniProt IDs to gene names
# region uniprot_to_gname
# uniprot = list(set(Reactome_BioGRID_additions['SubstrateID'].tolist()+Reactome_BioGRID_additions['KinaseID'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# uniprot = [x for x in uniprot if x not in uniprot_to_gname.keys()]
# for uni in uniprot:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+accession:{uni}&fields=gene_primary'
#     response = requests.get(query_url)
#     gene = response.text.split('\n')[1:-1]
#     uniprot_to_gname[uni] = gene
# with open(r'.\Variables\uniprot_to_gname.pkl', 'wb') as f:
#     pickle.dump(uniprot_to_gname, f)
# endregion
Reactome_BioGRID_additions['protein1'] = Reactome_BioGRID_additions.apply(lambda x: uniprot_to_gname[x['SubstrateID']][0] if uniprot_to_gname[x['SubstrateID']]!=[] else pd.NA, axis=1)
Reactome_BioGRID_additions['protein2'] = Reactome_BioGRID_additions.apply(lambda x: uniprot_to_gname[x['KinaseID']][0] if uniprot_to_gname[x['KinaseID']]!=[] else pd.NA, axis=1)
Reactome_BioGRID_additions.dropna(inplace=True)
Reactome_BioGRID_additions = Reactome_BioGRID_additions.iloc[:,[7,8]]
merge1 = pd.merge(breast_cancer_excl1, Reactome_BioGRID_additions, left_on=['source', 'target'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(breast_cancer_excl1, Reactome_BioGRID_additions, left_on=['source', 'target'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
breast_cancer_additions2 = pd.concat([merge1, merge2])
breast_cancer_additions2.drop_duplicates(inplace=True)
breast_cancer_additions2 = breast_cancer_additions2.iloc[:, :2]
breast_cancer_additions2.to_csv(r'..\Results\Network_Analysis\Reactome_BioGRID_additions(breast_cancer).tsv', sep="\t", index=False)

breast_cancer_prop_per = pd.DataFrame()
breast_cancer_prop_per = pd.concat([breast_cancer_prop_per,
                         pd.DataFrame({'Interactomes': ['Pathway Commons + BioGRID', 'BioGRID + Reactome'],
                                       '# Exclusive interactions': [len(breast_cancer_excl), len(breast_cancer_excl1)],
                                       '# Additional interactions (exclusive)': [len(breast_cancer_additions1), len(breast_cancer_additions2)],
                                       'Autophosphorylation': [len(breast_cancer_additions1[breast_cancer_additions1['source']==breast_cancer_additions1['target']]), len(breast_cancer_additions2[breast_cancer_additions2['source']==breast_cancer_additions2['target']])]}).set_index("Interactomes")])
# breast_cancer_prop_per.to_csv(r'..\Results\Network_Analysis\breast_cancer_prop_per.tsv', sep="\t")

## Cellphone DB
# Load the PPI file
cellphonedb = pd.read_csv(r"..\Data\v5.0.0_interaction_input.csv").iloc[:,1:] # Skip the first column
# Extract only PPIs
cellphonedb = cellphonedb[cellphonedb['is_ppi']==True]
# print(pd.concat([cellphonedb['partner_a'],cellphonedb['partner_b']]).nunique()) # 1140 proteins
cellphonedb = cellphonedb.iloc[:,:4]
cellphonedb.dropna(inplace=True)
# Get the gene names
cellphonedb['protein_name_a'] = cellphonedb['protein_name_a'].str.split('_HUMAN', expand=True)[0]
cellphonedb['protein_name_b'] = cellphonedb['protein_name_b'].str.split('_HUMAN', expand=True)[0]
#Map gene names to protein families to filter out those which have "kinase"
# region gname_to_pfam
# gene_names= list(set(cellphonedb['protein_name_a'].tolist()+cellphonedb['protein_name_b'].tolist()))
# base_url = 'https://rest.uniprot.org/uniprotkb'
# gene_names = [x for x in gene_names if x not in gname_to_pfam.keys()]
# for gene_name in gene_names:
#     query_url = f'{base_url}/search?format=tsv&query=organism_id:9606+gene:{gene_name}&fields=protein_families'
#     response = requests.get(query_url, timeout=10)
#     protein_name = response.text.split('\n')[1]
#     if "kinase" in protein_name.lower():
#         gname_to_pfam[gene_name] = protein_name
# with open('.\Variables\gname_to_pfam.pkl', 'wb') as f:
#     pickle.dump(gname_to_pfam, f)
# endregion
#Retain phosphorylation interactions
kinases = gname_to_pfam.keys()
cellphonedb_a = cellphonedb[cellphonedb['protein_name_a'].isin(kinases)]
cellphonedb_b = cellphonedb[cellphonedb['protein_name_b'].isin(kinases)]
cellphonedb = pd.concat([cellphonedb_a, cellphonedb_b])
cellphonedb.drop_duplicates(inplace=True)

# Check for interactions that are not in the chosen interactome networks
# PC + BioGRID
merge1 = pd.merge(cellphonedb, PC_new_BioGRID_excl, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(cellphonedb, PC_new_BioGRID_excl, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
cellphonedb_excl = pd.concat([merge1, merge2])
cellphonedb_excl.drop_duplicates(inplace=True)
cellphonedb_excl = cellphonedb_excl.iloc[:,:4]
# Additions that are contributed by CellphoneDB exclusively
merge1 = pd.merge(cellphonedb_excl, PC_BioGRID_additions, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(cellphonedb_excl, PC_BioGRID_additions, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
cellphonedb_additions1 = pd.concat([merge1, merge2])
cellphonedb_additions1.drop_duplicates(inplace=True)
cellphonedb_additions1 = cellphonedb_additions1.iloc[:, :4]
cellphonedb_additions1.to_csv(r'..\Results\Network_Analysis\PC_BioGRID_additions(CellphoneDB).tsv', sep="\t", index=False)

# BioGRID + Reactome
merge1 = pd.merge(cellphonedb, Reactome_BioGRID_excl, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(cellphonedb, Reactome_BioGRID_excl, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
cellphonedb_excl1 = pd.concat([merge1, merge2])
cellphonedb_excl1.drop_duplicates(inplace=True)
cellphonedb_excl1 = cellphonedb_excl1.iloc[:,:4]
# Additions that are contributed by CellphoneDB exclusively
merge1 = pd.merge(cellphonedb_excl1, Reactome_BioGRID_additions, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein1', 'protein2'], how='left', indicator=True)
merge1 = merge1[merge1['_merge'] == 'left_only'].drop(columns=['_merge'])
merge2 = pd.merge(cellphonedb_excl1, Reactome_BioGRID_additions, left_on=['protein_name_a', 'protein_name_b'], right_on=['protein2', 'protein1'], how='left', indicator=True)
merge2 = merge2[merge2['_merge'] == 'left_only'].drop(columns=['_merge'])
cellphonedb_additions2 = pd.concat([merge1, merge2])
cellphonedb_additions2.drop_duplicates(inplace=True)
cellphonedb_additions2 = cellphonedb_additions2.iloc[:, :4]
cellphonedb_additions2.to_csv(r'..\Results\Network_Analysis\Reactome_BioGRID_additions(CellphoneDB).tsv', sep="\t", index=False)

cellphonedb_prop_per = pd.DataFrame()
cellphonedb_prop_per = pd.concat([cellphonedb_prop_per,
                         pd.DataFrame({'Interactomes': ['Pathway Commons + BioGRID', 'BioGRID + Reactome'],
                                       '# Exclusive interactions': [len(cellphonedb_excl), len(cellphonedb_excl1)],
                                       '# Additional interactions (exclusive)': [len(cellphonedb_additions1), len(cellphonedb_additions2)],
                                       'Autophosphorylation': [len(cellphonedb_additions1[cellphonedb_additions1['protein_name_a']==cellphonedb_additions1['protein_name_b']]), len(cellphonedb_additions2[cellphonedb_additions2['protein_name_a']==cellphonedb_additions2['protein_name_b']])]}).set_index("Interactomes")])
# cellphonedb_prop_per.to_csv(r'..\Results\Network_Analysis\CellphoneDB_prop_per.tsv', sep="\t")

## Check if there are any overlap between the interactions in breast cancer cell lines and CellphoneDB
merge1 = pd.merge(breast_cancer, cellphonedb, how='inner', left_on=['source', 'target'], right_on=['protein_name_a', 'protein_name_b'])
merge2 = pd.merge(breast_cancer, cellphonedb, how='inner', left_on=['source', 'target'], right_on=['protein_name_b', 'protein_name_a'])
# There are no common interactions