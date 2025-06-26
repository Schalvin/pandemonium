import os
import subprocess
import argparse
from pathlib import Path
from statistics import mean
import networkx as nx
import pandas as pd
from natsort import natsorted
from sklearn_extra.cluster import KMedoids

def cluster_tb_stats(cluster_file, inflation):
    tm_scores = []
    identity = []
    length = []
    outliers = []
    nclusters = 0
    n_clusters_1 = 0
    with open(cluster_file) as file:
        for line in file:
            if line.startswith('cluster'):
                cluster = line.split(' : ')[-1].strip().split()
                if len(cluster) == 1:
                    n_clusters_1 += 1
                nclusters += 1
            if line.startswith('mean tm'):
                tm_scores.append(float(line.split(' : ')[-1].strip()))
            if line.startswith('mean identity'):
                identity.append(float(line.split(' : ')[-1].strip()))
            if line.startswith('mean length'):
                length.append(float(line.split(' : ')[-1].strip()))
            if line.startswith('outliers'):
                outliers.append(len((line.split(' : ')[-1].strip().split(' '))))
    mean_tm = sum(tm_scores)/(nclusters-n_clusters_1)
    mean_outliers = sum(outliers)/(nclusters-n_clusters_1)
    mean_identity = sum(identity)/(nclusters-n_clusters_1)
    mean_length = sum(length)/(nclusters-n_clusters_1)
    min_tm = min(tm_scores)
    min_outliers = min(outliers)
    min_identity = min(identity)
    min_length = min(length)
    max_tm = max(tm_scores)
    max_outliers = max(outliers)
    max_identity = max(identity)
    max_length = max(length)
    print(f'{inflation:<5s}{str(nclusters):<12s}{str(round(min_tm, 2)):<6s}{str(round(max_tm, 2)):<6s}{str(round(mean_tm, 2)):<15s}{str(round(min_identity, 2)):<6s}{str(round(max_identity, 2)):<6s}{str(round(mean_identity, 2)):<14s}{str(round(min_length, 2)):<8s}{str(round(max_length, 2)):<8s}{str(round(mean_length, 2)):<9s}{str(round(min_outliers, 2)):<6s}{str(round(max_outliers, 2)):<6s}{str(round(mean_outliers, 2)):<15s}{str(n_clusters_1):<24s}')


def clusters_table(cluster_files_format, input_dir, reverse = False):
    folder = Path(input_dir)
    cluster_files = folder.glob(f'*{cluster_files_format}*')
    if cluster_files_format.startswith('out'):
        cluster_files_inflation = [(str(cluster_file).split('.')[-1], cluster_file) for cluster_file in cluster_files]
        cluster_files_inflation = natsorted(cluster_files_inflation, key=lambda x: x[0], reverse = reverse)
        print('I    Nclusters')
    else:
        cluster_files_inflation = [(str(cluster_file).split('/')[-1].split('_')[0], cluster_file) for cluster_file in cluster_files]
        cluster_files_inflation = natsorted(cluster_files_inflation, key=lambda x: x[0], reverse = reverse)
        print('I    Nclusters   min   max   Avg TM-scores  min   max   Avg Identity  min     max     Avg Len  min   max   Avg outliers   Ncluster with one protein ')
    inflation_clusters_dict = {}
    for inflation, cluster_file in cluster_files_inflation:
        if cluster_files_format.startswith('out'):
            with open(cluster_file) as file:
                inflation_clusters_dict[inflation] = [cluster.split() for cluster in file]
                nclusters = len(inflation_clusters_dict[inflation])
                print(f'{inflation:<4s}, {str(nclusters):<9s}')
        else:
            cluster_tb_stats(cluster_file, inflation)
            with open(cluster_file) as file:
                inflation_clusters_dict[inflation] = [cluster.replace("cluster : ", "").split() for cluster in file if cluster.startswith("cluster : ")]
    print(inflation_clusters_dict)
    return inflation_clusters_dict

def create_graph(edge_file, node_file):
    node_table = pd.read_csv(node_file, sep='\t', usecols = ['SWORD2', 'length'], header=0)
    edge_table = pd.read_csv(edge_file, sep=' ', header=None, names=['Node1', 'Node2', 'T-score', 'Identity'])
    G = nx.Graph()
    # Add edges to the graph with the 'Value' as edge attribute
    for _, row in node_table.iterrows():
        G.add_node(row['SWORD2'], length=row['length'])
    for _, row in edge_table.iterrows():
        G.add_edge(row['Node1'], row['Node2'], tm_score=row['T-score'], identity=row['Identity'])
    return G

def get_mean_node_attribute_dict(graph, attribute):

    mean_values = {}
    n_nodes = len(graph.nodes())
    for node in graph.nodes():
        # Get the edges connected to the node
        edges = list(graph.edges(node, data=True))
        # If there are edges, calculate the mean of the edge weights
        if edges and n_nodes > 1:
            values = [data[attribute] for _, _, data in edges]
            mean_values[node] = sum(values) / (n_nodes - 1)
        else:
            mean_values[node] = 0
    # Output the mean value for all nodes 
    return mean_values

def get_mean_cluster_attribute_dict(graph, attribute):

    mean_values = {}
    nodes = graph.nodes()
    if len(nodes) > 1:
        print(nodes)
        values = [data[attribute] for _, data in nodes.items()]
        mean_value = sum(values) / (len(nodes) - 1)
    else:
        mean_value = 0
    return mean_value

def find_outliers(data_dict):
    import numpy as np

    values = np.array(list(data_dict.values()))
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1

    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    # Identify outliers
    outliers = [key for key, value in data_dict.items()
                if value < lower_bound or value > upper_bound]

    return outliers

def get_ref_structure_cluster(main_graph, cluster):
    #get dict with mean tm scores with all other proteins for each protein in the cluster
    cluster_graph = main_graph.subgraph(cluster)
    mean_t_dict = get_mean_node_attribute_dict(cluster_graph, 'tm_score')
    mean_identity_dict = get_mean_node_attribute_dict(cluster_graph, 'identity')
    mean_length = get_mean_cluster_attribute_dict(cluster_graph, 'length')
    # get outlier values :
    outliers = find_outliers(mean_t_dict)
    #calculate 3 kmedoids for given cluster 
    medoids = []
    if len(cluster) > 3:
        adj_matrix = nx.adjacency_matrix(cluster_graph, nodelist = cluster).toarray()
        medoid_index = KMedoids(n_clusters=3).fit(adj_matrix).medoid_indices_
        for index in medoid_index:
            medoids.append(cluster[index])
    return max(mean_t_dict, key=mean_t_dict.get), medoids, outliers, mean(mean_t_dict.values()), mean(mean_identity_dict.values()), mean_length

def relate_ref_structures_between_clusterings(cluster_files_format, input_dir, main_graph_file, nodes_files, inflations = None, output='tree.graphml'):
    # this function uses a graph (main_graph), where all nodes (prot domain) are connected and have a TM score value, 
    # and a list of nodes from a cluster file (mcl) to find the reference of each cluster based on the highest mean tm score for the cluster
    main_graph = create_graph(main_graph_file, nodes_files)

    adj_mat = pd.DataFrame(nx.adjacency_matrix(main_graph, weight='tm_score').toarray())
    adj_mat.index = main_graph.nodes()
    adj_mat.columns = main_graph.nodes()

    print(f'number of proteins : {main_graph.number_of_nodes()}')
    inflation_clusters_dict = clusters_table(cluster_files_format, input_dir, reverse = True)
    inflation_clusters_ref_dict = {}
    if not inflations :
        inflations = [inflation for inflation, clusters in inflation_clusters_dict.items()]
    else: inflations = inflations.split(',')
    inflations = natsorted(inflations, reverse = True)
    print(inflations[-1])
    #data for UMAP :
    adj_mat.to_csv(f"umap_adjmatrix_{'_'.join(inflations)}.tsv", sep = "\t")
    # mapper = umap.UMAP().fit(adj_mat)
    # uplot.points(mapper)
    for inflation, clusters in inflation_clusters_dict.items():
        cluster_ref_list = [(cluster, get_ref_structure_cluster(main_graph, cluster)) for cluster in clusters]
        inflation_clusters_ref_dict[inflation]=cluster_ref_list

    #output all kmedoids
    for inflation, clusters in inflation_clusters_ref_dict.items():
        with open(f"{inflation}_cluster_kmedoids_outliers.txt", 'w') as kmedoids_outfile:
            for cluster, (reference, kmedoids, outliers, mean_tm, mean_identity, mean_length) in clusters:
                kmedoids_outfile.write(f'reference : {reference}\n')
                kmedoids_outfile.write(f'mean tm : {mean_tm}\n')
                kmedoids_outfile.write(f'mean length : {mean_length}\n')
                kmedoids_outfile.write(f'mean identity : {mean_identity}\n')
                kmedoids_outfile.write(f'kmedoids : {" ".join(kmedoids)}\n')
                kmedoids_outfile.write(f'cluster : {" ".join(cluster)}\n')
                if len(outliers) != 0:
                    kmedoids_outfile.write(f'outliers : {" ".join(outliers)}\n')

    #create hierarchical tree with cluster references
    G = nx.DiGraph()

    for i in range(len(inflations)-1):
        print(inflations[i+1], inflations[i])
        for cluster_1, (reference_1, _, _, _, _, _) in inflation_clusters_ref_dict[inflations[i]]:
            for cluster_2, (reference_2, _, _, _, _, _) in inflation_clusters_ref_dict[inflations[i+1]]:
                if reference_1 in cluster_2:
                    G.add_node(f'{reference_1}_{inflations[i]}', size=len(cluster_1))
                    if f'{reference_2}_{inflations[i+1]}' not in G.nodes():
                        G.add_node(f'{reference_2}_{inflations[i+1]}', size=len(cluster_2))
                    G.add_edge(f'{reference_2}_{inflations[i+1]}', f'{reference_1}_{inflations[i]}')
    
    G.add_node('artificial_central_node', size=10)
    print('artificial', inflations[-1], len(inflation_clusters_ref_dict[inflations[-1]]))
    for cluster, (reference, _, _, _, _, _) in inflation_clusters_ref_dict[inflations[-1]]:
        if f'{reference}_{inflations[-1]}' not in G.nodes():
            G.add_node(f'{reference}_{inflations[-1]}', size=len(cluster))
        G.add_edge('artificial_central_node', f'{reference}_{inflations[-1]}')
    nx.write_graphml(G, output)

def align_queries_to_target(target, queries, folder, max_queries):
    print(target)
    log_file = os.path.join(target)
    queries = [os.path.join(folder, str(x).strip() + '.pdb') for x in queries]
    if max_queries and len(queries) > max_queries:
        queries = queries[:max_queries]
    target = os.path.join(folder, target) + '.pdb'
    if len(queries) > 1 :
        command = ["kpax", "-multi", "-nopivot", "-unified"] + [target] + queries
    else: 
        command = ["kpax", "-unified"] + [target] + queries
        
    with open(log_file, 'w') as log:
        try:
            subprocess.run(command, check=True, stdout=log)
        except subprocess.CalledProcessError as e:
            print(f"Error while executing command: {e}")

def align_hierarchical_clustering(tree_file, input_dir, output_dir, max_queries = None):
    tree = nx.read_graphml(tree_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)
    for target in tree.nodes():
        queries = [query for query in tree.successors(target)]
        if queries:
            align_queries_to_target(target, queries, input_dir, max_queries)
    return

def align_kmedoids(kmedoids_file, input_dir, output_dir, align_outliers = False):
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)
    with open(kmedoids_file) as kfile:
        for line in kfile.readlines():
            if line.startswith('reference'):
                reference = line.split(' : ')[-1].strip()
            elif line.startswith('kmedoids'):
                kmedoids = line.split(' : ')[-1].strip()
                if kmedoids:
                    kmedoids = kmedoids.split()
                    align_queries_to_target(reference, kmedoids, input_dir, max_queries = 3)
            elif line.startswith('outliers') and align_outliers:
                kmedoids = line.split(' : ')[-1].strip()
                kmedoids = kmedoids.split()
                os.chdir(os.path.join(output_dir, 'outliers'))
                align_queries_to_target(reference, outliers)
                os.chdir(output_dir)
                continue 

def prot_cluster_tb(prot_tb, cluster_file_format, inflations):
    prot_tb = pd.read_csv(prot_tb, sep='\t', header=0)



def combine_prot_tb_clusters(input_dir, prot_tb_file, cluster_file_format, inflations = None):
    prot_tb = pd.read_csv(prot_tb_file, sep='\t', header=0)
    inflation_clusters_dict = clusters_table(cluster_file_format, input_dir, reverse = True)
    print(inflation_clusters_dict.keys())
    if not inflations :
        inflations = [inflation for inflation, clusters in inflation_clusters_dict.items()]
    else: inflations = inflations.split(',')

    for inflation in inflations:
        for idx, cluster in enumerate(inflation_clusters_dict[inflation]):
            for prot in cluster:
                prot_tb.loc[prot_tb['SWORD2'] == prot, inflation] = str(idx)
    prot_tb = prot_tb.dropna()
    prot_tb.to_csv(f'{prot_tb_file.replace(".tsv", "")}_mcmlclusters.tsv', sep = "\t", index = None)
    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--programm", type=str, help="'table', 'graph_references', 'align_references', 'align_kmedoids'")
    parser.add_argument("-c", "--cluster_files", type=str, help="cluster files name format if multiple cluster files, otherwise path to the cluster file")
    parser.add_argument("-i", "--input_dir", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-n", "--node_file", type=str)
    parser.add_argument("-e", "--edge_file", type=str)
    parser.add_argument("-I", "--inflations", type=str, help="inflation values to keep (format : I17 for 1.5, I50 for 5), comma separated")
    parser.add_argument("-t", "--tree", type=str)
    parser.add_argument("-k", "--kmedoids_file", type=str)
    parser.add_argument("-m", "--max_queries", type=int, default=None)
    parser.add_argument("--prot_tb", type=str)


    args = parser.parse_args()
    if args.programm == "table":
        clusters_table(args.cluster_files, args.input_dir)
    if args.programm == "graph_references":
        relate_ref_structures_between_clusterings(args.cluster_files, args.input_dir, args.edge_file, args.node_file, args.inflations, args.output)
    if args.programm == "align_references":
        align_hierarchical_clustering(args.tree, args.input_dir, args.output, args.max_queries)
    if args.programm == "align_kmedoids":
        align_kmedoids(args.kmedoids_file, args.input_dir, args.output)
    if args.programm == "align_outliers":
        align_kmedoids(args.kmedoids_file, args.input_dir, args.output, align_outliers = True)
    if args.programm == "prot_cluster_tb":
        combine_prot_tb_clusters(args.input_dir, args.prot_tb, args.cluster_files, args.inflations)

