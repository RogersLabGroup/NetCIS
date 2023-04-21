#!/usr/bin/python
# Mathew Fischbach, fisch872@umn.edu

from pathlib import Path
from multiprocessing import Pool

from docopt import docopt
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns

import utils


def load_args():
    doc = """
    Usage: main.py --dir DIR [options]

     --dir=DIR                      path to data
     --meta=FILE                    tab delimited file for providing additional information per input file per line

    Options:
     -h --help                      show this help message and exit
     --verbose=N                    print more verbose information using 1, 2, or 3 [default: 0]
     --output=DIR                   a directory to save results to, else results are saved to default directory [default: output]
     --jobs=N                       number of processes to run [default: 1]
     --threshold=N                  maximum distance to connect two insertions together in the CIS network [default: 50000]
     --rand_state=N                 seed for random states [default: 42]
     --loglevel=LEVEL               set logging level for information and plots. Can us debug, info, warning, error, critical [default: warning]
    """

    input_args = docopt(doc)

    new_input_args = {}
    for key, value in input_args.items():
        new_key = key.split('-')[-1]
        new_input_args[new_key] = value

    int_opts = ['jobs', 'rand_state', 'verbose', 'threshold']
    for opts in int_opts:
        new_input_args[opts] = int(new_input_args[opts])

    new_input_args['loglevel'] = new_input_args['loglevel'].upper()

    new_input_args['dir'] = Path(new_input_args['dir'])
    new_input_args['meta'] = Path(new_input_args['meta'])
    new_input_args['output'] = Path(new_input_args['output'])

    return new_input_args

def create_graph(insertions_df, chrom, threshold):
    G = nx.Graph()
    chrom_df = insertions_df[insertions_df['chr'] == chrom].sort_values('pos')
    
    # for each insertion
    for i in range(len(chrom_df.index)):
        insert = chrom_df.iloc[i]
        new_node = insert['pos']
        # add node(i) as an insertion location into the network
        G.add_node(new_node,
                counts=insert['count'],
                counts_irr=insert['count_irr'],
                counts_irl=insert['count_irl'],
                orient=insert['tpn promoter orient'],
                chrom=insert['chr'],
                position=insert['pos'],
        )
        
        for other_node in G.nodes:
            # find distance between nodes using their position
            node_dist = abs(other_node - new_node)
            # don't add edge to self
            if node_dist == 0:
                continue
            # if distance between node(i) and node(j) is less than threshold
            if node_dist <= threshold:
                # add edge(ij) with a weight of the distance (or inverse?) to the network
                G.add_edge(new_node, other_node, weight=1/node_dist)

    return G

def graph_properties(G):
    print(f"number of nodes: {G.number_of_nodes()}")
    print(f"number of edges: {G.number_of_edges()}")
    num_inserts = 0
    for node in G.nodes:
        num_inserts += G.nodes[node]['counts']
    print(f"number of insertions: {num_inserts}")

def main(args):
    # logger = utils.get_logger(args['output'], log_level=args['loglevel'])

    ### Load data
    # meta info about each file
    meta_df = pd.read_csv(args['meta'])

    # get all files in data dir and load each file as pandas.DataFrame
    insert_list = []
    for file in args['dir'].iterdir():
        file_meta_info = meta_df.loc[meta_df['file name'] == file.name].drop('file name', axis=1)
        tmp_df = pd.read_csv(file)
        tmp_df['cell type'] = file_meta_info['cell type'].tolist()[0]
        tmp_df['cell id'] = file_meta_info['cell id'].tolist()[0]
        tmp_df['tumor type'] = file_meta_info['tumor type'].tolist()[0]
        insert_list.append(tmp_df)
    inserts_df = pd.concat(insert_list, ignore_index=True)

    # separate data into case/controls and combine insertions
    # Group at the chr, pos, and tpn promoter orientation level. Keep count, count_irr and count_irl
    group_cols = ['chr', 'pos', 'tpn promoter orient']
    keep_cols = ['count', 'count_irr', 'count_irl']

    insert_case = inserts_df[inserts_df['tumor type'] == 'S']
    insert_case_grouped = insert_case.groupby(by=group_cols, as_index=False, dropna=False).sum(numeric_only=True)[group_cols + keep_cols]

    insert_control = inserts_df[inserts_df['tumor type'] != 'S']
    insert_control_grouped = insert_control.groupby(by=group_cols, as_index=False, dropna=False).sum(numeric_only=True)[group_cols + keep_cols]

    ### Construct network (using pseudo code from graph framework article)
    # TODO: this is where we could parallelize
    # for each node in the network
    # with Pool(njobs) as p:
    # iter_args = [ (test_node, insert['pos'], threshold, G) for test_node in G.nodes ]
    # p.starmap(create_graph_helper, iter_args)
    # p.starmap_async(create_graph_helper, iter_args)
    for chrom in np.unique(inserts_df['chr'].to_numpy()):
        chr1_case_graph = create_graph(insert_case_grouped, chrom, args['threshold'], args['njobs'])
        nx.write_graphml(chr1_case_graph, args['output'] / f'case-{chrom}.graphml')

        chr1_control_graph = create_graph(insert_control_grouped, chrom, args['threshold'], args['njobs'])
        nx.write_graphml(chr1_control_graph, args['output'] / f'control-{chrom}.graphml')

    ### Analyze network
    #     export network into CytoScape and analyze there OR use networkx for analysis
    #     for each connected subgraph in the graph
    #         if subgraph is not a random network
    #             add subgraph to list of non-random CIS
    #     Explore non-random CIS
    #     Somehow get a p-value of each CIS


if __name__ == '__main__':
    main(load_args())
