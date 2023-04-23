from pathlib import Path
from multiprocessing import Pool
from typing import Generator

from docopt import docopt
import numpy as np
from pandas import read_csv, concat, DataFrame
import networkx as nx


def load_args() -> dict:
    doc = """
    Generate common insertion sites (CIS) using a network (graph) based approach. 
    The only parameter that will change how a CIS is generated can be control with --threshold
    
    Usage: cis_networks.py --output_prefix DIR [options]
    
     --output_prefix=DIR            a prefix of the output directory that will have "-graphs" appended to it

    Options:
     -h --help                      show this help message and exit
     --verbose=N                    print more verbose information using 1 or 2 [default: 0]
     --jobs=N                       number of processes to run [default: 1]
     --threshold=N                  maximum distance to connect two insertions together in the CIS network. We suggest not going over the default value [default: 50000]
    """

    input_args = docopt(doc)

    new_input_args = {}
    for key, value in input_args.items():
        new_key = key.split("-")[-1]
        new_input_args[new_key] = value

    int_opts = ["jobs", "verbose", "threshold"]
    for opts in int_opts:
        new_input_args[opts] = int(new_input_args[opts])

    new_input_args["insertion_dir"] = Path(new_input_args["output_prefix"] + "-insertions")
    new_input_args["output"] = Path(new_input_args["output_prefix"] + "-graphs")

    return new_input_args

def create_graph(chrom_df: DataFrame, threshold, save_file, verbose=0) -> None:
    G = nx.Graph()
    
    # prepare the insertions by grouping them together
    # find the total count of insertions and the counts per sequencing library (IRR/IRL)
    insert_cols = ['chr', 'pos', 'tpn promoter orient', 'seq library']
    tmp = chrom_df.groupby(by=insert_cols, as_index=False, dropna=False)['read name'].count()
    tmp['count'] = tmp.pop('read name')
    count_irr = np.where(tmp['seq library'] == 'IRR', tmp['count'], 0)
    count_irl = np.where(tmp['seq library'] == 'IRL', tmp['count'], 0)
    tmp.insert(5, "count_irr", count_irr)
    tmp.insert(6, "count_irl", count_irl)
    
    # group insertions without the sequencing library. 
    # As long as the transposon orientation, chromosome, and position are the same, 
    # then it does not matter which library the insertion came from
    node_cols = ['chr', 'pos', 'tpn promoter orient']
    insertion_nodes = tmp.groupby(by=node_cols, as_index=False, dropna=False).sum(numeric_only=True)
    insertion_nodes['read names'] = chrom_df.groupby(by=node_cols, dropna=False, group_keys=False)['read name'].apply(list).reset_index(drop=True)
    
    # TODO: for some reason there are few insertions that occur both in IRR and IRL
    # both_libs = insertion_nodes[ (insertion_nodes['count_irl'] != 0) & (insertion_nodes['count_irr'] != 0) ]
    
    # add each insertion. Since they are unique, I can add the edges after all the nodes are in
    for i in range(len(insertion_nodes)):
        if (i % 1000 == 0) and (i != 0) and verbose:
            print(f"\t{i+1/len(insertion_nodes)} insertions")
        insert = insertion_nodes.iloc[i]
        new_node = f"{insert['pos']}|{insert['tpn promoter orient']}"
        # add node(i) as an insertion location into the network
        G.add_node(
            new_node,
            counts=insert["count"],
            counts_irr=insert["count_irr"],
            counts_irl=insert["count_irl"],
            orient=insert["tpn promoter orient"],
            chrom=insert["chr"],
            position=insert["pos"],
        )
        # for other_node in G.nodes:
        #     if other_node == new_node:
        #         continue
        #     # find distance between nodes using their position
        #     node_dist = abs(G.nodes[other_node]["position"] - G.nodes[new_node]["position"])
        #     # double check don't add edge to self
        #     if node_dist == 0:
        #         continue
        #     # if distance between node(i) and node(j) is less than threshold
        #     if node_dist <= threshold:
        #         # add edge(ij) with a weight of the distance (or inverse?) to the network
        #         
        #         G.add_edge(new_node, other_node, weight=1 / node_dist)
    
    # the following code does what is commented out above but I am keeping all of this in for a future user to reference
    
    # nodes are inherently ordered as they are added in the graph. 
    # however, the ordering doens't have to numerically make sense
    ordered_nodes = G.nodes()
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
    cis_nodes = dist_nodes < threshold  # symmetric 2d array
    
    # get the indices of the lower left triangle of the symmetric matrix.
    # edges_ind is a tuple of two array. The same index location in both arrays is used 
    # to index a single value from the symmetric matrix. This results in two very long 
    # arrays that will index all the values of the lower left triangle of the matrix
    edges_ind = np.tril_indices_from(cis_nodes, k=-1) # tuple of two 1d arrays
    
    # keep nodes that are under the threshold
    keep_nodes = cis_nodes[edges_ind]  # 1d array
    
    # set up the nodes to be a numpy array for easy indexing
    ordered_nodes = np.array(G.nodes())  # 1d array
    
    # get the actual node names for the lower left triangle via as the column
    nodes1 = ordered_nodes[edges_ind[1][keep_nodes]]  # 1d array
    # the rows
    nodes2 = ordered_nodes[edges_ind[0][keep_nodes]]  # 1d array
    # and edge weights (TODO: which can be modified for a differnt weighting method, maybe 1 / log10(x) instead?)
    nodes_dist = 1 / dist_nodes[edges_ind][keep_nodes]  # 1d array
    # combine the nodes and weights into an iterable that can be passed wholly into the graph
    # an edge is defined as the first node, the second node, and then a dict of attributes, such as weight
    edges_to_add = [ (x, y, {"weight": z}) for x, y, z in zip(nodes1, nodes2, nodes_dist) ]
    G.add_edges_from(edges_to_add)

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
        cell_type, cell_id, tumor_type = file.stem.split("-")
        tmp_df = read_csv(file)
        tmp_df["cell type"] = cell_type
        tmp_df["cell id"] = cell_id
        tmp_df["tumor type"] = tumor_type
        insert_list.append(tmp_df)
    inserts_df = concat(insert_list, ignore_index=True)

    # separate data into case/controls
    insert_case = inserts_df[inserts_df["tumor type"] == "S"]
    insert_control = inserts_df[inserts_df["tumor type"] != "S"]

    # get all chromosomes to separate further the case/controls dataframes
    chrom_list = np.unique(inserts_df["chr"].to_numpy())
    # don't allow more jobs than there are chromosomes
    jobs = args["jobs"]
    if len(chrom_list) < jobs:
        print(f"Reducing number of jobs from {jobs} to {len(chrom_list)}, since there are only {len(chrom_list)} chromosomes present.")
        jobs = len(chrom_list)
        
    # construct CIS network per chromosome for case and control insertions
    iter_gen = create_graph_generator(chrom_list, insert_case, insert_control, args['threshold'], out_dir_case, out_dir_control)
    with Pool(jobs) as p:
        [ x for x in p.imap_unordered(create_graph_helper, iter_gen) ]
        
if __name__ == "__main__":
    main(load_args())
