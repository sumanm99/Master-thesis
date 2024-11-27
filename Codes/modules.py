'''
Some user-defined functions
'''
import pandas as pd
import ast

def find_overlaps(t1, t2):
    '''
    Find the overlaps between UniProt IDs that are present as tuple of strings

    Input: pandas.core.series.Series object
    t1: column of dataframe 1 (tuple of strings)
    t2: column of dataframe 2 (tuple of strings)

    Return UniProt ID overlaps as tuple of strings
    '''
    #Convert pandas.core.series.Series to tuple of strings
    t1 = ast.literal_eval(t1) if isinstance(t1, str) else tuple(t1)
    t2 = ast.literal_eval(t2) if isinstance(t2, str) else tuple(t2)
    return tuple(set(t1).intersection(set(t2)))

def find_overlaps2(s1, t2):
    '''
    Find the overlaps between UniProt IDs that are present in tuple of strings and string

    Input: pandas.core.series.Series object
    s1: column of dataframe 1 (string)
    t2: column of dataframe 2 (tuple of strings)

    Return UniProt ID overlaps as string
    '''
    return s1 if any(x in s1 for x in t2) else None

def check_condition(df1, df2):
    '''
    Checks if the Substrate ID and Kinase ID (strings) of dataframe df2 are present in 'protein1' and 'protein2' (tuple of strings) of dataframe df1, respectively.

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    matches = df2.loc[(df2['SubstrateID'].isin(df1['protein1'])) & (df2['KinaseID'].isin(df1['protein2']))]
    if matches.empty:
        return pd.DataFrame()
    else:
        return pd.concat([df1.to_frame().T]*len(matches)).reset_index(drop=True).join(matches.reset_index(drop=True))

def check_condition_reverse(df1, df2):
    '''
    Checks if the Substrate ID and Kinase ID (strings) of dataframe df2 are present in 'protein2' and 'protein1' (tuple of strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    matches = df2.loc[(df2['SubstrateID'].isin(df1['protein2'])) & (df2['KinaseID'].isin(df1['protein1']))]
    if matches.empty:
        return pd.DataFrame()
    else:
        return pd.concat([df1.to_frame().T]*len(matches)).reset_index(drop=True).join(matches.reset_index(drop=True))

def check_condition2(df1, df2):
    '''
    Checks if the Substrate ID (string) and Kinase ID (tuple of strings) of dataframe df2 overlap with 'protein1' and 'protein2' (strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    # merge df1 and df2 on common columns
    merged_df = df1.merge(df2, left_on='protein1', right_on='SubstrateID')
    # filter out rows where protein2 is not in KinaseID
    merged_df = merged_df[merged_df.apply(lambda x: x['protein2'] in x['KinaseID'], axis=1)]
    return merged_df

# #Checking the function
# df1 = pd.DataFrame({
#     'protein1': ['A', 'B', 'C'],
#     'protein2': ['Z', 'X', 'Z']
# })
# df2 = pd.DataFrame({
#     'SubstrateID': ['C', 'D', 'E'],
#     'KinaseID': [('Z', 'X'), ('W', 'Y'), ('X', 'Z')]
# })
# merge = check_condition2(df1, df2)

def check_condition2_reverse(df1, df2):
    '''
    Checks if the Substrate ID (string) and Kinase ID (tuple of strings) of dataframe df2 overlap with 'protein2' and 'protein1' (strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    # merge df1 and df2 on common columns
    merged_df = df1.merge(df2, left_on='protein2', right_on='SubstrateID')
    # filter out rows where protein1 is not in KinaseID
    merged_df = merged_df[merged_df.apply(lambda x: x['protein1'] in x['KinaseID'], axis=1)]
    return merged_df

def check_condition3(df1, df2):
    '''
    Checks if the Substrate ID (string) and Kinase ID (tuple of strings) of dataframe df2 are present in 'protein1' (string) and 'protein2' (tuple of strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    # merge df1 and df2 on common columns
    merged_df = df1.merge(df2, left_on='protein1', right_on='SubstrateID')
    # filter out rows where protein2 and KinaseID do not have any overlaps
    merged_df = merged_df[merged_df.apply(lambda x: any(item in x['protein2'] for item in x['KinaseID']), axis=1)]
    return merged_df

def check_condition3_reverse(df1, df2):
    '''
    Checks if the Substrate ID (string) and Kinase ID (tuple of strings) of dataframe df2 are present in 'protein2' (string) and 'protein1' (tuple of strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    # merge df1 and df2 on common columns
    merged_df = df1.merge(df2, left_on='protein2', right_on='SubstrateID')
    # filter out rows where protein1 and KinaseID do not have any overlaps
    merged_df = merged_df[merged_df.apply(lambda x: any(item in x['protein1'] for item in x['KinaseID']), axis=1)]
    return merged_df

def check_condition4(df1, df2):
    '''
    Checks if the Substrate ID (string) and Kinase ID (tuple of strings) of dataframe df2 are present in 'protein1' and 'protein2' (tuple of strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    merged_data = []
    for _, row1 in df1.iterrows():
        for _, row2 in df2.iterrows():
            if any(item in row2['SubstrateID'] for item in row1['protein1']) and any(item in row2['KinaseID'] for item in row1['protein2']):
                merged_data.append(list(row1) + list(row2))
    merged_df = pd.DataFrame(merged_data, columns=df1.columns.tolist() + df2.columns.tolist())
    return merged_df

# #Checking the function
# df1 = pd.DataFrame({
#     'protein1': [('C', 'A'), ('W', 'Y')],
#     'protein2': [('Z', 'X'), ('W', 'Y')]
# })
# df2 = pd.DataFrame({
#     'SubstrateID': ['C', 'D', 'E'],
#     'KinaseID': [('Z', 'X'), ('W', 'Y'), ('X', 'Z')]
# })
# merge = check_condition4(df1, df2)

def check_condition4_reverse(df1, df2):
    '''
    Checks if the Substrate ID (string) and Kinase ID (tuple of strings) of dataframe df2 are present in 'protein2 and 'protein1' (tuple of strings) of dataframe df1, respectively

    Inputs:
    df1: dataframe that contains 'protein1' and 'protein2' columns (Interactomes)
    df2: dataframe that contains 'SubstrateID' and 'KinaseID' columns

    Return the overlaps as a dataframe
    '''
    merged_data = []
    for _, row1 in df1.iterrows():
        for _, row2 in df2.iterrows():
            if any(item in row2['SubstrateID'] for item in row1['protein1']) and any(item in row2['KinaseID'] for item in row1['protein2']):
                merged_data.append(list(row1) + list(row2))
    merged_df = pd.DataFrame(merged_data, columns=df1.columns.tolist() + df2.columns.tolist())
    return merged_df