from pathlib import Path
from multiprocessing import Pool
from typing import Generator
import sys

from tqdm import tqdm
from docopt import docopt
import numpy as np
from pandas import read_csv, concat, DataFrame
import networkx as nx

from analysis import graph_properties

    
def load_args() -> dict:
    doc = """  
    Generate common insertion sites (CIS) using a network (graph) based approach. 
    The only parameter that will change how a CIS is generated can be control with --threshold
    
    Usage: cis_networks.py --output_prefix DIR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-graphs" appended to it

    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   (TODO: how to allow --verbose meaning 1 as well as supplying value?) print more verbose information using 0, 1 or 2 [default: 0]
     -t, --threshold=N                 maximum distance to connect two insertions together in the CIS network. We suggest not going over the default value [default: 50000]
     -j, --jobs=N                      number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "threshold", "jobs"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])

    new_args["insertion_dir"] = Path(new_args["output_prefix"] + "-insertions")
    new_args["output"] = Path(new_args["output_prefix"] + "-graphs")

    return new_args

def add_nodes(insertion_df):
    def add_node(insert):
        node = insert.pos
        attr = {
            "position": insert.pos,
            "chrom": insert.chr,
            "counts": insert.counts,
            "counts_irr": insert.counts_irr,
            "counts_irl": insert.counts_irl,
            "counts_trp_orient_pos": insert.counts_trp_orient_pos,
            "counts_trp_orient_neg": insert.counts_trp_orient_neg,
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
        # double check don't add edge to self
        if node_dist == 0:
            continue
        # if distance between node(i) and node(j) is less than threshold
        if node_dist <= threshold:
            # add edge(ij) with a weight of the distance (or inverse?) to the network
            
            G.add_edge(new_node, other_node, weight=1 / node_dist)
    """ 
    # nodes are inherently ordered as they are added in the graph,
    # however, the ordering doens't have to numerically make sense for this function
    
    
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

def create_graph(chrom_df: DataFrame, save_file, threshold=50000, verbose=0) -> None:
    G = nx.Graph()
    chrom_df.insert(4, "counts_irr", np.where(chrom_df['library'] == 'IRR', 1, 0))
    chrom_df.insert(5, "counts_irl", np.where(chrom_df['library'] == 'IRL', 1, 0))
    chrom_df.insert(6, "counts_trp_orient_pos", np.where(chrom_df['tpn_promoter_orient'] == '+', 1, 0))
    chrom_df.insert(7, "counts_trp_orient_neg", np.where(chrom_df['tpn_promoter_orient'] == '-', 1, 0))
    cols = ["counts_irr", "counts_irl", "counts_trp_orient_pos", "counts_trp_orient_neg"]
    
    tmp_group = chrom_df.groupby(by=['chr', 'pos'], sort=False, as_index=False, dropna=False)
    insertion_nodes_df = tmp_group[cols].sum()
    insertion_nodes_df.insert(2, "counts", tmp_group['read_name'].count().pop('read_name'))

    # add nodes and edges to graph
    G.add_nodes_from(add_nodes(insertion_nodes_df))
    G.add_edges_from(find_edges(G.nodes(), threshold))
    
    if verbose > 1:
        graph_properties(G)

    # save the graph
    nx.write_graphml(G, save_file)

def create_graph_helper(iter_args) -> None:
    insert_case_chrom, case_file, insert_control_chrom, control_file, threshold, verbose = iter_args
    create_graph(insert_case_chrom, case_file, threshold, verbose)
    create_graph(insert_control_chrom, control_file, threshold, verbose)
    
def create_graph_generator(chrom_list, insert_case, insert_control, case_dir, control_dir, args) -> Generator[tuple, None, None]:
    for chrom in chrom_list:
        insert_case_chrom = insert_case[insert_case['chr'] == chrom]    
        insert_control_chrom = insert_control[insert_control['chr'] == chrom]
        case_file = case_dir / f"{chrom}.graphml"
        control_file = control_dir / f"{chrom}.graphml"
        yield ( insert_case_chrom, case_file, insert_control_chrom, control_file, args["threshold"], args["verbose"] )

def main(args) -> None:
    # prepare output
    out_dir_case = args['output'] / "case"
    out_dir_case.mkdir(parents=True, exist_ok=True)
    out_dir_control = args['output'] / "control"
    out_dir_control.mkdir(parents=True, exist_ok=True)
    
    # get all files in data dir, load each file as pandas.DataFrame, and add meta data based on the file name
    insert_list = []
    for file in args["insertion_dir"].iterdir():
        tmp_df = read_csv(file, sep="\t")
        tumor_model, sample_id, tissue_type = file.name.split("-")
        tmp_df["tumor_model"] = tumor_model
        tmp_df["sample_id"] = sample_id
        tmp_df["tissue"] = tissue_type.split(".")[0]  # RT/LT/S
        insert_list.append(tmp_df)
    inserts_df = concat(insert_list, ignore_index=True)


    # TODO: how are we choosing case and controls?
    # TODO: create another script to compare the union vs intersection of insertions between left and right tumors
    # also then look at and compare insertions counts of intersection of insertions and the not intersection (Laura, 5/11/23)
    # or more formally, the disjointed insertions
    
    
    # TODO: average the left and right insertions, if they are at the same insertion site
    # separate data into case/controls
    insert_case = inserts_df[inserts_df["tissue"] != "S"]
    insert_control = inserts_df[inserts_df["tissue"] == "S"]

    # get all chromosomes to separate further the case/controls dataframes
    chrom_list = np.unique(inserts_df["chr"].to_numpy())
    
    # don't allow more jobs than there are chromosomes
    jobs = args["jobs"]
    num_chr = len(chrom_list)
    if num_chr < jobs:
        print(f"Reducing number of jobs from {jobs} to {num_chr}, since there are only {num_chr} chromosomes present.")
        jobs = len(chrom_list)
        
    # construct CIS network per chromosome for case and control insertions
    iter_gen = create_graph_generator(chrom_list, insert_case, insert_control, out_dir_case, out_dir_control, args)
    iter_gen = tqdm(iter_gen, total=num_chr)
    with Pool(jobs) as p:
        for _ in p.imap_unordered(create_graph_helper, iter_gen):
            pass
        p.close()
        
if __name__ == "__main__": 
    main(load_args())
    