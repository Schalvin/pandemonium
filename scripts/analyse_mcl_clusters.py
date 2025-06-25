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
    mean_tm = sum(tm_scores)/(nclusters-n_clusters_1)
    print(f'{inflation:<4s}{str(nclusters):<9s}{str(round(mean_tm, 2)):<14s}{str(n_clusters_1):<25s}')


def clusters_table(cluster_files_format, input_dir, reverse = False):
    folder = Path(input_dir)
    cluster_files = folder.glob(f'*{cluster_files_format}*')
    if cluster_files_format.startswith('out'):
        cluster_files_inflation = [(str(cluster_file).split('.')[-1], cluster_file) for cluster_file in cluster_files]
        cluster_files_inflation = natsorted(cluster_files_inflation, key=lambda x: x[0], reverse = reverse)
        print('I   Nclusters')
    else:
        cluster_files_inflation = [(str(cluster_file).split('/')[-1].split('_')[0], cluster_file) for cluster_file in cluster_files]
        cluster_files_inflation = natsorted(cluster_files_inflation, key=lambda x: x[0], reverse = reverse)
        print('I   Nclusters  Mean TM-scores  Ncluster with one protein')
    inflation_clusters_dict = {}
    for inflation, cluster_file in cluster_files_inflation:
        if cluster_files_format.startswith('out'):
            with open(cluster_file) as file:
                inflation_clusters_dict[inflation] = [cluster.split() for cluster in file]
                nclusters = len(inflation_clusters_dict[inflation])
                print(f'{inflation:<4s}, {str(nclusters):<9s}')
        else:
            cluster_tb_stats(cluster_file, inflation)
    return inflation_clusters_dict

def create_graph(graph_file):
    graph_table = pd.read_csv(graph_file, sep=' ', header=None, names=['Node1', 'Node2', 'T-score'])
    G = nx.Graph()
    # Add edges to the graph with the 'Value' as edge attribute
    for _, row in graph_table.iterrows():
        G.add_edge(row['Node1'], row['Node2'], weight=row['T-score'])
    return G

def get_mean_node_tm_dict(graph):

    mean_values = {}
    n_nodes = len(graph.nodes())
    for node in graph.nodes():
        # Get the edges connected to the node
        edges = list(graph.edges(node, data=True))
        # If there are edges, calculate the mean of the edge weights
        if edges and n_nodes > 1:
            values = [data['weight'] for _, _, data in edges]
            mean_values[node] = sum(values) / (n_nodes - 1)
        else:
            mean_values[node] = 0
    # Output the mean value for all nodes 
    return mean_values

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
    mean_t_dict = get_mean_node_tm_dict(cluster_graph)
    # get outlier values :
    outliers = find_outliers(mean_t_dict)
    #calculate 3 kmedoids for given cluster 
    medoids = []
    if len(cluster) > 3:
        adj_matrix = nx.adjacency_matrix(cluster_graph, nodelist = cluster).toarray()
        medoid_index = KMedoids(n_clusters=3).fit(adj_matrix).medoid_indices_
        for index in medoid_index:
            medoids.append(cluster[index])
    return max(mean_t_dict, key=mean_t_dict.get), medoids, outliers, mean(mean_t_dict.values())

def relate_ref_structures_between_clusterings(cluster_files_format, input_dir, main_graph_file, inflations = None, output='tree.graphml'):
    # this function uses a graph (main_graph), where all nodes (prot domain) are connected and have a TM score value, 
    # and a list of nodes from a cluster file (mcl) to find the reference of each cluster based on the highest mean tm score for the cluster
    main_graph = create_graph(main_graph_file)
    print(f'number of proteins : {main_graph.number_of_nodes()}')
    inflation_clusters_dict = clusters_table(cluster_files_format, input_dir, reverse = True)
    inflation_clusters_ref_dict = {}
    if not inflations :
        inflations = [inflation for inflation, clusters in inflation_clusters_dict.items()]
    else: inflations = inflations.split(',')
    inflations = natsorted(inflations, reverse = True)
    print(inflations[-1])
    for inflation, clusters in inflation_clusters_dict.items():
        cluster_ref_list = [(cluster, get_ref_structure_cluster(main_graph, cluster)) for cluster in clusters]
        inflation_clusters_ref_dict[inflation]=cluster_ref_list

    #output all kmedoids
    for inflation, clusters in inflation_clusters_ref_dict.items():
        with open(f"{inflation}_cluster_kmedoids_outliers.txt", 'w') as kmedoids_outfile:
            for cluster, (reference, kmedoids, outliers, mean_tm) in clusters:
                kmedoids_outfile.write(f'reference : {reference}\n')
                kmedoids_outfile.write(f'mean tm : {mean_tm}\n')
                kmedoids_outfile.write(f'kmedoids : {" ".join(kmedoids)}\n')
                kmedoids_outfile.write(f'cluster : {" ".join(cluster)}\n')
                if len(outliers) != 0:
                    kmedoids_outfile.write(f'outliers : {" ".join(outliers)}\n')

    #create hierarchical tree with cluster references
    G = nx.DiGraph()

    for i in range(len(inflations)-1):
        print(inflations[i+1], inflations[i])
        for cluster_1, (reference_1, _, _, _) in inflation_clusters_ref_dict[inflations[i]]:
            for cluster_2, (reference_2, _, _, _) in inflation_clusters_ref_dict[inflations[i+1]]:
                if reference_1 in cluster_2:
                    G.add_node(f'{reference_1}_{inflations[i]}', size=len(cluster_1))
                    if f'{reference_2}_{inflations[i+1]}' not in G.nodes():
                        G.add_node(f'{reference_2}_{inflations[i+1]}', size=len(cluster_2))
                    G.add_edge(f'{reference_2}_{inflations[i+1]}', f'{reference_1}_{inflations[i]}')
    
    G.add_node('artificial_central_node', size=10)
    print('artificial', inflations[-1], len(inflation_clusters_ref_dict[inflations[-1]]))
    for cluster, (reference, _, _, _) in inflation_clusters_ref_dict[inflations[-1]]:
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--programm", type=str, help="'table', 'graph_references', 'align_references', 'align_kmedoids'")
    parser.add_argument("-c", "--cluster_files", type=str, help="cluster files name format if multiple cluster files, otherwise path to the cluster file")
    parser.add_argument("-i", "--input_dir", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-g", "--graph", type=str)
    parser.add_argument("-I", "--inflations", type=str, help="inflation values to keep (format : I17 for 1.5, I50 for 5), comma separated")
    parser.add_argument("-t", "--tree", type=str)
    parser.add_argument("-k", "--kmedoids_file", type=str)
    parser.add_argument("-m", "--max_queries", type=int, default=None)

    args = parser.parse_args()
    if args.programm == "table":
        clusters_table(args.cluster_files, args.input_dir)
    if args.programm == "graph_references":
        relate_ref_structures_between_clusterings(args.cluster_files, args.input_dir, args.graph, args.inflations, args.output)
    if args.programm == "align_references":
        align_hierarchical_clustering(args.tree, args.input_dir, args.output, args.max_queries)
    if args.programm == "align_kmedoids":
        align_kmedoids(args.kmedoids_file, args.input_dir, args.output)
    if args.programm == "align_outliers":
        align_kmedoids(args.kmedoids_file, args.input_dir, args.output, align_outliers = True)

