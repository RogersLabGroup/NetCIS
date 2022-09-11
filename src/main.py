#!/usr/bin/python
# Mathew Fischbach, fisch872@umn.edu

from pathlib import Path

from docopt import docopt
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns

from . import utils


def load_args():
    doc = """
    Usage: main.py --dir DIR [options]
    
     --dir=DIR                   path to data
     
    Options:
     -h --help                   show this help message and exit
     --verbose=N                 print more verbose information using 1, 2, or 3 [default: 0]
     --output=DIR                a directory to save results to, else results are saved to default directory [default: output]
     --jobs=N                    number of processes to run [default: 1]
     --rand_state=N              seed for random states [default: 42]
     --loglevel=LEVEL            set logging level for information and plots. Can us debug, info, warning, error, critical [default: warning]
    """
    
    input_args = docopt(doc)
    
    new_input_args = {}
    for key, value in input_args.items():
        new_key = key.split('-')[-1]
        new_input_args[new_key] = value
    
    int_opts = ['jobs', 'rand_state', 'verbose']
    for opts in int_opts:
        new_input_args[opts] = int(new_input_args[opts])
    
    new_input_args['loglevel'] = new_input_args['loglevel'].upper()
    
    return new_input_args


def main():
    args = load_args()
    output_path = Path(args['output'])
    log_level = args['loglevel']
    logger = utils.get_logger(output_path, log_level)
    
    
    ### Load data
    data_dir = Path(args['dir'])
    # get all files in data dir
    
    # load each file
    
    # convert to numpy data structures


    ### Construct network (using pseudo code from graph framework article)
    #     Split insertion sites (IS) by chromosome
    #     Order IS for each chromosome
    #     for each IS(i)
    #         add node(i) as an IS location into the network
    #         for each node (j) in the network
    #             if distance between node(i) and node(j) is less than threshold
    #                 add edge(ij) with a weight of the distance (or inverse?) to the network  
    #     split transcriptional element (TE) by chromosomes
    #     order TE for each chromosome 
    #     for each TE(i)
    #         add node(i) as a TE location into the network
    #         for each node (j) in the graph
    #             if distance between node(i) and node(j) is less than threshold
    #                 add edge(ij) with a weight of the distance (or inverse?) to the network
        
        
    ### Analyze network 
    #     export network into CytoScape and analyze there OR use networkx for analysis
    #     for each connected subgraph in the graph
    #         if subgraph is not a random network
    #             add subgraph to list of non-random CIS
    #     Explore non-random CIS
    #     Somehow get a p-value of each CIS


if __name__ == '__main__':
    
    # Preprocess data before running this script
    #     cutadapt to trim adapter and SB tags
    #     map reads to reference with bowtie2
    #     output as bam files
    #     *change bam files to bed files, one insertion site per line?
    
    main()
    
