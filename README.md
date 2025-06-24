## Protein domain classification using SWORD2 domains and clustering based on Kpax alignemnts

### Steps followed :
#### Identification of similar domains in the ATLAS protein database (~2000 proteins)
- Step 1 : cutting SWORD2 domains out from protein pdb files and output domain correspondance table with ECOD, CATH and SCOP domains (/scripts/cut_pdb.py)
- Step 2 : Align pair by pair each ATLAS domain to the others with Kpax (scripts/kpax_allpair-align.py) 
- Step 3 : From Kpax log files, extract all TM scores to directed weighted graph with query, target, TM-score, filter edges to get undirected graph with only edges with TM-scores > 0.45 and query and target len > 40 (keep only one edge for a pair of proteins, the one for which the TM-score is highest) and output an both in abc format (scripts/kpax_to_T_graph.py)
- Step 4 : use MCL clustering on the abc adjacency table, with TM-score as similarity, using different inflation values (scripts/mcl_clusterings_ATLAS.sh)

#### Choice of reference domains in ATLAS protein database
- Step 5 : Analyse mcl clusters, define reference structures, 3 medoids, outliers for each cluster, output a directed graph relating all cluster references (scripts/analyse_mcl_clusters.py)