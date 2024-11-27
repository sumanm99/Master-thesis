## Parsing_the_input_files.py  
- **Drop rows with NAs to only retain hits with complete information in the standard format (UniProt ID of kinases were missing).** 
  - In GPS, 1836 UniProt kinase IDs were missing out of the 16626 entries (edges/rows) for Homo sapiens.  
  - In NetworKIN, all kinases had Ensembl protein IDs, so there was no need to drop NAs. The conversion to UniProt IDs was done for only 4038 out of 7143 entries.  
  - All of the GPS peptides are 61-mers, but not all peptides in the NetworKIN file are 9-mers (marked in red in NetworKIN Excel file).  
- Finding overlaps,  
  - The slight variations in the positions in the predictions of GPS and NetworKIN could be due to the fact that they were trained on datasets from different databases.  
  - Substrate, Kinase, Code and Peptide (s1) - The position variations were ignored. Positions represent the correct phosphorylation sites. But the position numbers are slightly different in the available predictions.  
  - **From s1, I noticed that one of the hits had a large difference in the position numbers, so I limited the difference to 2.** 
  - The motif and context score of NetworKIN were included, and sorted in descending order to view the most confident predictions. The (min, max) of motif score and content score were (0.472, 1) and (0.139, 0.999), respectively. **The average value of motif and context scores were around 0.576 and 0.876, respectively, which seem to be around the threshold of 0.5 and 0.9 as mentioned in the paper.  
  - **All the result files have the motif scores above the threshold of 0.5.** I filtered out the predictions with context scores <0.9. Very few predictions were lost.  
  - **I could not find confidence scores for the predictions made by GPS (I've mailed the author).**

Check if the phosphorylated proteins in prostate cancer is predicted by these tools.  
- The curation of the Prostate cancer based predictions (Stathmin) in the Excel sheet was based on the following,  
  - For NetworKIN, I searched for the keyword "Stathmin" in the "description_substrate" column.  
  - For GPS, I curated the corresponding Substrate IDs referring to the NetworKIN curated list.  
  - There are no hits for stathmin in the overlaps of the prediction tools, that are of high confidence.    
  - **P13668 and P21818 belong to Rat, but it has been listed under Homo sapiens in GPS prediction file. I excluded them after manually verifying the organism with UniProt.**  
  - **O70166 belongs to Mouse, but it has been presented in NetworKIN prediction file. The stathmin UniProt list does not contain this ID. I noticed this after comparing the results with the Prostate cancer excel sheets.** 
  - There are fewer hits in NetworKIN because of the formatted file being used; not all EnsemblID has a corresponding UniProt ID.  
Check for overlaps between the predictions and the databases that they have used for relevance.  
- NetworKIN and Phospho.ELM database
  - The NetworKIN prediction file is the human phosphorylation network (HPN).  
  - The Phospho.ELM curation was obtained from [link](http://phospho.elm.eu.org/dataset.html) for vertebrates.  
  - We have about 37145 kinase-substrate interactions for *Homo sapiens* in the Phospho.ELM file.    
  - Some of the substrates are represented by Ensembl protein IDs. 277 predictions were dropped which contained the unconverted Ensembl protein IDs.  
  - The entries with missing kinases were not dropped because it is not used to find the overlaps.  
  - **The kinases in Phospho.ELM don't seem to have any pattern with 'predictor' or 'genesymbol_kinase' of NetworKIN, to check for correspondence.**  
- GPS and PhosphoSitePlus database  
  - The curation file from PhosphoSitePlus database was obtained from [link](https://www.phosphosite.org/staticDownloads). "Kinase_Substrate_Dataset.gz" file has been used in the GPS paper as well.  
  - There are 13828 kinase-substrate interactions for *Homo sapiens* in the PhosphoSitePlus file. Only the entries with both the kinase organism and the substrate organism coming from humans were retained.  

Venn diagrams  
- The diagram is constructed based on the number of predictions (edges).
- The vertices and edges are represented as (vertices, edges).
- **For the Phospho.ELM and NetworKIN diagram, the number of vertices of the overlap is not accurate, but only approximate, due to incomparable kinase info.**
