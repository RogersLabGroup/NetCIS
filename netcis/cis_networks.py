from pathlib import Path
from multiprocessing import Pool
from typing import Generator
import sys

from tqdm import tqdm
from docopt import docopt
import numpy as np
from pandas import read_csv, concat, DataFrame
import networkx as nx

    
def load_args() -> dict:
    doc = """  
    Generate common insertion sites (CIS) using a network (graph) based approach. 
    The only parameter that will change how a CIS is generated can be control with --threshold
    
    Usage: cis_networks.py --output_prefix DIR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-graphs" appended to it

    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information using 0, 1 or 2 [default: 0]
     -t, --threshold=N                 maximum distance to connect two insertions together in the CIS network. We suggest not going over the default value [default: 50000]
     -j, --jobs=N                      number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { new_args[key.split("-")[-1]]: value for key, value in docopt(doc).items() }

    int_opts = ["jobs", "verbose", "threshold"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])

    new_args["insertion_dir"] = Path(new_args["output_prefix"] + "-insertions")
    new_args["output"] = Path(new_args["output_prefix"] + "-graphs")

    return new_args

def add_nodes(insertion_df):
    def add_node(insert):
        node = f"{insert.pos}|{insert.tpn_promoter_orient}"
        attr = {
            "counts": insert.count,
            "counts_irr": insert.count_irr,
            "counts_irl": insert.count_irl,
            "orient": insert.tpn_promoter_orient,
            "chrom": insert.chr,
            "position": insert.pos,
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


    # remove the transposon orientation from the end of the node name
    tmp_order = [ int(x.split("|")[0]) for x in ordered_nodes ]
    # check if this changes the number of unique nodes.
    # If we have + and - at the same location, this assert will fail.
    # This isn't a bad thing but I want to know when it is happening
    assert len(np.unique(ordered_nodes)) == len(np.unique(tmp_order))
    
    # cast the nodes into a numpy array that can be used to broadcast into a symmetric matrix of distances
    nodes = np.array(tmp_order).reshape(-1, 1)
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
    
    # set up the nodes to be a numpy array for easy indexing
    ordered_nodes = np.array(ordered_nodes)  # 1d array
    
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

def create_graph(chrom_df: DataFrame, threshold, save_file, verbose=0) -> None:
    G = nx.Graph()
    
    # prepare the insertions by grouping them together
    # find the total count of insertions and the counts per sequencing library (IRR/IRL)
    insert_cols = ['chr', 'pos', 'tpn_promoter_orient', 'library']
    tmp = chrom_df.groupby(by=insert_cols, as_index=False, dropna=False)['read_name'].count()
    tmp['count'] = tmp.pop('read name')
    count_irr = np.where(tmp['library'] == 'IRR', tmp['count'], 0)
    count_irl = np.where(tmp['library'] == 'IRL', tmp['count'], 0)
    tmp.insert(5, "count_irr", count_irr)
    tmp.insert(6, "count_irl", count_irl)
    
    # group insertions without the sequencing library. 
    # As long as the transposon orientation, chromosome, and position are the same, 
    # then it does not matter which library the insertion came from
    node_cols = ['chr', 'pos', 'tpn_promoter_orient']
    insertion_nodes = tmp.groupby(by=node_cols, as_index=False, dropna=False).sum(numeric_only=True)
    insertion_nodes['read_names'] = chrom_df.groupby(by=node_cols, dropna=False, group_keys=False)['read_name'].apply(list).reset_index(drop=True)
    
    # TODO: for some reason there are few insertions that occur both in IRR and IRL. 
    # Why is that and does this change with the new preprocessing scripts?
    # both_libs = insertion_nodes[ (insertion_nodes['count_irl'] != 0) & (insertion_nodes['count_irr'] != 0) ]
    
    # add nodes and edges to graph
    G.add_nodes_from(add_nodes(insertion_nodes))
    G.add_edges_from(find_edges(G.nodes(), threshold))

    # save the graph
    nx.write_graphml(G, save_file)

def create_graph_helper(iter_args) -> None:
    insert_case_chrom, insert_control_chrom, threshold, case_file, control_file = iter_args
    create_graph(insert_case_chrom, threshold, case_file)
    create_graph(insert_control_chrom, threshold, control_file)
    
def create_graph_generator(chrom_list, insert_case, insert_control, threshold, case_dir, control_dir) -> Generator[tuple, None, None]:
    for chrom in chrom_list:
        print(chrom)
        insert_case_chrom = insert_case[insert_case['chr'] == chrom]    
        insert_control_chrom = insert_control[insert_control['chr'] == chrom]
        case_file = case_dir / f"{chrom}.graphml"
        control_file = control_dir / f"{chrom}.graphml"
        yield ( insert_case_chrom, insert_control_chrom, threshold, case_file, control_file )

def main(args) -> None:
    # prepare output
    out_dir_case = args['output'] / "case"
    out_dir_case.mkdir(parents=True, exist_ok=True)
    out_dir_control = args['output'] / "control"
    out_dir_control.mkdir(parents=True, exist_ok=True)
    
    # get all files in data dir, load each file as pandas.DataFrame, and add meta data based on the file name
    insert_list = []
    for file in args["insertion_dir"].iterdir():
        tmp_df = read_csv(file)  # TODO: load .tsv
        
        tumor_model, sample_id, tumor_tmp = file.name.split("-")
        tissue_type, lib_tmp = tumor_tmp.split("_")
        library = lib_tmp.split(".")[0]
    
        tmp_df["tumor_model"] = tumor_model
        tmp_df["sample_id"] = sample_id
        tmp_df["tissue"] = tissue_type  # RT/LT/S
        tmp_df["library"] = library  # IRR/IRL
    
        insert_list.append(tmp_df)
        
    inserts_df = concat(insert_list, ignore_index=True)

    # TODO: how are we choosing case and controls?
    # TODO: create another script to compare the union vs intersection of insertions between left and right tumors
    # also then look at and compare insertions counts of intersection of insertions and the not intersection (Laura, 5/11/23)
    # or more formally, the disjointed insertions
    
    
    # separate data into case/controls
    insert_case = inserts_df[inserts_df["tissue"] == "S"]
    insert_control = inserts_df[inserts_df["tissue"] != "S"]

    # get all chromosomes to separate further the case/controls dataframes
    chrom_list = np.unique(inserts_df["chr"].to_numpy())
    # don't allow more jobs than there are chromosomes
    jobs = args["jobs"]
    if len(chrom_list) < jobs:
        print(f"Reducing number of jobs from {jobs} to {len(chrom_list)}, since there are only {len(chrom_list)} chromosomes present.")
        jobs = len(chrom_list)
        
    # construct CIS network per chromosome for case and control insertions
    iter_gen = create_graph_generator(chrom_list, insert_case, insert_control, args['threshold'], out_dir_case, out_dir_control)
    iter_gen = tqdm(iter_gen)
    with Pool(jobs) as p:
        [ x for x in p.imap_unordered(create_graph_helper, iter_gen) ]
        p.close()
        
if __name__ == "__main__": 
    main(load_args())
    