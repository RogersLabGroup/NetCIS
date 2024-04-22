from pathlib import Path
from multiprocessing import Pool
import pickle

from tqdm import tqdm
from docopt import docopt
import numpy as np
import pandas as pd
from pandas import read_csv, concat, DataFrame
import networkx as nx

    
def load_args() -> dict:
    doc = """  
    Generate pseudo common insertion sites (pCIS) using a network (graph) based approach. 
    The only parameter that will change how a CIS is generated can be controled with --threshold
    
    Usage: 
        pcis_networks.py --output_prefix DIR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-graphs" appended to it

    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information if available using 0, 1 or 2 [default: 0]
     -t, --threshold=N                 maximum distance to connect two insertions together in the CIS network. We suggest not going over the default value [default: 50000]
     -j, --njobs=N                     number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "threshold", "njobs"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])

    new_args["insertion_dir"] = Path(new_args["output_prefix"] + "-insertions")
    new_args["depth_dir"] = Path(new_args["output_prefix"] + "-insertions-depth")
    new_args["output"] = Path(new_args["output_prefix"] + "-graphs")
    
    if new_args["verbose"] > 1:
        print(new_args)
        
    return new_args

def graph_properties(G, verbose=0):
    """
    Calculate various properties of a graph.

    Args:
        G (networkx.Graph): The input graph
        verbose (int, optional): Verbosity level. Default is 0.

    Returns:
        dict: Dictionary containing graph properties
    """
    nodes = G.number_of_nodes()
    edges = G.number_of_edges()
    num_inserts = sum([ G.nodes[node]["counts"] for node in G.nodes ])
    subgraphs_by_nodes = sorted(nx.connected_components(G), key=len, reverse=True)
    num_subgraphs = len(subgraphs_by_nodes)
    if verbose:
        print(f"number of nodes: {nodes}")
        print(f"number of edges: {edges}")
        print(f"number of insertions: {num_inserts}")
        print(f"number of subgraphs (pCIS) {num_subgraphs}")
    return {"nodes": nodes, "edges": edges, "num_inserts": num_inserts, "num_subgraphs": num_subgraphs}

def add_nodes(insertion_df):
    def add_node(insert):
        node = insert.pos
        attr = {
            "position": insert.pos,
            "chrom": insert.chr,
            "CPM": insert.CPM,
            "counts": insert.counts,
            "counts_irr": insert.counts_irr,
            "counts_irl": insert.counts_irl,
            # "counts_trp_orient_pos": insert.counts_trp_orient_pos,
            # "counts_trp_orient_neg": insert.counts_trp_orient_neg,
            "sample_IDs": insert.sample_IDs,
        }
        return (node, attr)
    return [ add_node(x) for x in insertion_df.itertuples() ]

def find_edges(ordered_nodes, threshold):
    """
    use of numpy vector based methods to speed up the process of identifying and adding edges 
    to a graph when using distance between nodes (for which the nodes are numbers).
    
    A notebook with the step-by-step output of this function is available under notebooks/create_graph_edges.ipynb
    
    This function does the same this, and is a direct replacement, for the following code 

    # find and add edges
    for other_node in G.nodes:
        if other_node == new_node:
            continue
        # find distance between nodes using their position
        node_dist = abs(G.nodes[other_node]["position"] - G.nodes[new_node]["position"])
        # double check don"t add edge to self
        if node_dist == 0:
            continue
        # if distance between node(i) and node(j) is less than threshold
        if node_dist <= threshold:
            # add edge(ij) with a weight of the distance (or inverse?) to the network
            
            G.add_edge(new_node, other_node, weight=1 / node_dist)
    """ 
    # nodes are inherently ordered as they are added in the graph,
    # however, the ordering doens"t have to numerically make sense for this function
    
    
    # set up the nodes to be a numpy array for easy indexing
    ordered_nodes = np.array(ordered_nodes)  # 1d array
    
    # reshape nodes to be used to broadcast into a symmetric matrix of distances
    nodes = ordered_nodes.reshape(-1, 1)
    dist_nodes = np.abs(nodes - nodes.T)  # symmetric 2d array
    
    # cis nodes are those that are under the threshold
    cis_nodes = dist_nodes <= threshold  # symmetric 2d array
    
    # get the indices of the lower left triangle of the symmetric matrix.
    # edges_ind is a tuple of two array. The same index location in both arrays is used 
    # to index a single value from the symmetric matrix. This results in two very long 
    # arrays that will index all the values of the lower left triangle of the matrix
    edges_ind = np.tril_indices_from(cis_nodes, k=-1) # tuple of two 1d arrays
    
    # keep nodes that are under the threshold
    keep_nodes = cis_nodes[edges_ind]  # 1d array
    
    # get the actual node names for the lower left triangle via as the column
    nodes1 = ordered_nodes[edges_ind[1][keep_nodes]]  # 1d array
    # the rows
    nodes2 = ordered_nodes[edges_ind[0][keep_nodes]]  # 1d array
    # and edge weights (TODO: which can be modified for a differnt weighting method, maybe 1 / log10(x) instead?)
    nodes_dist = 1 / dist_nodes[edges_ind][keep_nodes]  # 1d array
    # combine the nodes and weights into an iterable that can be passed wholly into the graph
    # an edge is defined as the first node, the second node, and then a dict of attributes, such as weight
    edges_to_add = [ (x, y, {"weight": z}) for x, y, z in zip(nodes1, nodes2, nodes_dist) ]
    return edges_to_add

def create_graph(iter_args):
    chrom_df, save_dir, threshold, verbose = iter_args
    
    G = nx.Graph()
    cols = ["CPM", "counts_irr", "counts_irl"]
    tmp_group = chrom_df.groupby(by=['chr', 'pos'], sort=False, as_index=False, dropna=False)
    insertion_nodes_df = tmp_group[cols].sum()
    insertion_nodes_df.insert(2, "counts", tmp_group['count'].count().pop('count'))

    # add in info about which samples are in each insertion site
    tmp_samples = chrom_df.groupby(by=['chr', 'pos'], sort=False, as_index=False, dropna=False)["sampleID"].apply(lambda x: x.unique())
    if tmp_samples.size == 0:
        insertion_nodes_df["n_samples"] = 0
        insertion_nodes_df["sample_IDs"] = []
    else:
        insertion_nodes_df.insert(6, "n_samples", tmp_samples["sampleID"].apply(lambda x: len(x)))
        insertion_nodes_df.insert(6, "sample_IDs", tmp_samples["sampleID"].apply(lambda x: list(x)).to_list())

    # add nodes and edges to graph
    G.add_nodes_from(add_nodes(insertion_nodes_df))
    G.add_edges_from(find_edges(G.nodes(), threshold))
    
    if verbose > 1:
        graph_properties(G)

    # save the graph
    nx.write_gml(G, save_dir / "G.gml")
    
    # save subgraphs from graph
    subgraphs_by_nodes = sorted(nx.connected_components(G), key=len, reverse=True)
    subgraphs = [ G.subgraph(x) for x in subgraphs_by_nodes ]
    with open(save_dir / "subgraphs.pickle", "wb") as f:
        pickle.dump(subgraphs, f, pickle.HIGHEST_PROTOCOL)
    
def create_graph_generator(chrom_list, treatment_inserts, treatment_dir, args):
    for chrom in chrom_list:
        treatment_chrom_inserts = treatment_inserts[treatment_inserts["chr"] == chrom]
        treatment_chrom_dir = treatment_dir / f"{chrom}"
        treatment_chrom_dir.mkdir(parents=True, exist_ok=True)
        
        yield ( treatment_chrom_inserts, treatment_chrom_dir, args["threshold"], args["verbose"] )

def main(args) -> None:
    # get all files in data dir, load each file as pandas.DataFrame
    insertion_list = [ pd.read_csv(file, sep="\t") for file in args["depth_dir"].iterdir() ]
    inserts_df = pd.concat(insertion_list, ignore_index=True)
    inserts_df.insert(4, "counts_irr", np.where(inserts_df['library'] == 'IRR', 1, 0))
    inserts_df.insert(5, "counts_irl", np.where(inserts_df['library'] == 'IRL', 1, 0))
    
    chrom_list = np.unique(inserts_df["chr"].to_numpy())
    treatment_list = inserts_df["treatment"].unique()

    # total unique samples across all treatments
    total_samples = inserts_df["sampleID"].unique().shape[0]
    metadata = {"total": total_samples}

    for treatment in treatment_list:
        print(treatment)
        # prepare output
        out_dir = args['output'] / treatment
        out_dir.mkdir(parents=True, exist_ok=True)
        
        treatment_df = inserts_df[inserts_df["treatment"] == treatment]
        metadata[treatment] = treatment_df["sampleID"].unique().shape[0]
        
        # don't allow more jobs than there are chromosomes
        jobs = args["njobs"]
        num_chr = len(chrom_list)
        if num_chr < jobs:
            print(f"Reducing number of jobs from {jobs} to {num_chr}, since there are only {num_chr} chromosomes present.")
            jobs = len(chrom_list)
            
        # construct CIS network per chromosome for treatment insertion
        iter_gen = create_graph_generator(chrom_list, treatment_df, out_dir, args)
        iter_gen = tqdm(iter_gen)
        with Pool(jobs) as p:
            for _ in p.imap_unordered(create_graph, iter_gen):
                pass
            p.close()
            
    # save sample numbers as meta data for network analysis
    samples, counts = zip(*metadata.items())
    meta_df = pd.DataFrame({"samples": samples, "counts": counts})
    meta_df.to_csv(args['output'].parent / "samples_with_insertions.csv", index=False)

if __name__ == "__main__": 
    main(load_args())
    