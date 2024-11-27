import pandas as pd
from unipressed import IdMappingClient
import time

#Load saved variables/results
overlaps = pd.read_csv('..\Results\overlaps_phosphorylation.tsv', sep='\t')
GPS_format = pd.read_csv('.\Variables\GPS_format.tsv', sep='\t')
NetworKIN_format = pd.read_csv(r'.\Variables\NetworKIN_format.tsv', sep='\t')

##Prostate cancer related predictions
#Stathmin
request = IdMappingClient.submit(
    source="Gene_Name", dest="UniProtKB", ids={"STMN1", "STMN2"}
)
time.sleep(5.0)
# print(list(request.each_result()))
stathmin_uniprot = list(map(lambda x: x['to'], list(request.each_result())))
mask = overlaps['SubstrateID'].isin(stathmin_uniprot)
stathmin = overlaps[mask]
#No hits

mask = GPS_format['SubstrateID'].isin(stathmin_uniprot)
stathmin_GPS = GPS_format[mask]
#There are fewer hits because of the formatted file being used.
mask = NetworKIN_format['SubstrateID'].isin(stathmin_uniprot)
stathmin_NetworKIN = NetworKIN_format[mask]

#All of the stathmin based predictions in standard format
stathmin_all = pd.concat([stathmin_GPS,stathmin_NetworKIN])
order = [0,3,2,1,4,5,6,7,8]
stathmin_all = stathmin_all[[stathmin_all.columns[i] for i in order]]
stathmin_all.sort_values(by=['SubstrateID', 'Position'], inplace=True)
stathmin_all.sort_values(by=['motif_score', 'context_score'], inplace=True, ascending=False)
stathmin_all = stathmin_all[(stathmin_all['context_score']>=0.9) | (stathmin_all['context_score'].isna())] #Filter out the predictions with context scores <0.9
stathmin_all.to_csv(r'..\Results\all_stathmin_hits.tsv', sep="\t", index=False)