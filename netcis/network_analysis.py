import pickle
from pathlib import Path
from multiprocessing import Pool

from docopt import docopt
import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import binomtest, ranksums, fisher_exact, boschloo_exact
from tqdm import tqdm


def load_args() -> dict:
    """
    Load command-line arguments using docopt.

    Returns:
        dict: Parsed command-line arguments
    """
    
    doc = """  
    Generate common insertion sites (CIS) using a network (graph) based approach. 
    The only parameter that will change how a CIS is generated can be control with --threshold

    Usage: 
        analysis.py --output_prefix DIR --ta_dir DIR --gene_annot FILE --case STR --control STR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-analysis" appended to it
     -e, --ta_dir=DIR                  directory that contains the TA locations for each chromosome in bed format
     -g, --gene_annot=FILE             MGI's mouse menetic markers excluding withdrawn genes
     -a, --case=STR                    treatment type value to use as case
     -b, --control=STR                 treatment type value to use as control
     
    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information using 0, 1 or 2 [default: 0]
     -t, --ta_error=N                  how many bases to expand the search for a TA site at each insertion [default: 5]
     -p, --pval_threshold=N            p-value to exclude pCIS for significance [default: 0.05]
     -x, --gene_expander=N             number of base pairs to extend the boundaries of genes when annotation genes to CIS [default: 50000]
     -j, --njobs=N                     number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "ta_error", "njobs", "gene_expander"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])
    
    new_args["pval_threshold"] = float(new_args["pval_threshold"])
    
    new_args["graph_dir"] = Path(new_args["output_prefix"] + "-graphs")
    new_args["ta_dir"] = Path(new_args["ta_dir"])
    new_args["gene_annot"] = Path(new_args["gene_annot"])
    new_args["output"] = Path(new_args["output_prefix"] + "-analysis")
    new_args["output"].mkdir(exist_ok=True)

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
    norm_num_inserts = sum([ G.nodes[node]['CPM'] for node in G.nodes ])
    subgraphs_by_nodes = sorted(nx.connected_components(G), key=len, reverse=True)
    num_subgraphs = len(subgraphs_by_nodes)
    if verbose:
        print(f"number of nodes: {nodes}")
        print(f"number of edges: {edges}")
        print(f"number of normliazed insertions: {norm_num_inserts}")
        print(f"number of subgraphs (pCIS) {num_subgraphs}")
    return {"nodes": nodes, "edges": edges, "norm_num_inserts": norm_num_inserts, "num_subgraphs": num_subgraphs}
    
def subgraph_properties(G, verbose=0):
    """
    Calculate properties of a subgraph.

    Args:
        G (networkx.Graph): The input subgraph
        verbose (int, optional): Verbosity level. Default is 0.

    Returns:
        dict: Dictionary containing subgraph properties
    """
    nodes = G.number_of_nodes()
    edges = G.number_of_edges()
    norm_num_inserts = sum([ G.nodes[node]['CPM'] for node in G.nodes ])
    tmp_pos = sorted([ G.nodes[node]["position"] for node in G.nodes ])
    min_pos = min(tmp_pos)
    max_pos = max(tmp_pos)
    range_pos = max_pos - min_pos
    
    sample_IDs = { x for y in [ G.nodes[node]["sample_IDs"] for node in G.nodes ] for x in y }
    num_unique_samples = len(sample_IDs)
    
    if verbose:
        print(f"number of nodes: {nodes}")
        print(f"number of edges: {edges}")
        print(f"number of normalized insertions: {norm_num_inserts}")
        print(f"min position: {min_pos}")
        print(f"max position: {max_pos}")
        print(f"range: {range_pos}")
        print(f"# of unique samples: {num_unique_samples}")
        
    return {"nodes": nodes, 
            "edges": edges, 
            "norm_num_inserts": norm_num_inserts, 
            "min_pos": min_pos, 
            "max_pos": max_pos, 
            "range": range_pos,
            "sample_IDs": [sample_IDs],
            "num_unique_samples": num_unique_samples,
            }
    
def subgraph_TA_sites(G, bed, ta_error, verbose=0):
    """
    Calculate properties related to TA sites in a subgraph.

    Args:
        G (networkx.Graph): The input subgraph
        bed (pandas.DataFrame): Bed file data containing TA site information
        ta_error (int): The error margin for TA site matching
        verbose (int, optional): Verbosity level. Default is 0.

    Returns:
        dict: Dictionary containing TA site properties
    """
    num_insert_sites = G.number_of_nodes()
    
    tmp_pos = sorted([ G.nodes[node]["position"] for node in G.nodes ])
    ta_sites = bed[(bed[1] > min(tmp_pos)) & (bed[2] < max(tmp_pos))]
    num_ta_sites = len(ta_sites)
    
    arr1 = np.array(tmp_pos).reshape(-1, 1)
    arr2 = ta_sites[1].to_numpy().reshape(-1, 1)
    arr3 = ta_sites[1].to_numpy().reshape(-1, 1)
    ta_inserts = (arr1 >= (arr2.T - ta_error)) & (arr1 <= (arr3.T + ta_error))
    num_ta_insert_sites = ta_inserts.any(axis=1).sum()
    
    if verbose:
        print(f"number of insertions in subgraph: {num_insert_sites}")
        print(f"number of TA sites in subgraph: {num_ta_sites}")
        print(f"number of insertions within a TA site (+/- {ta_error} bp): {num_ta_insert_sites}")
    
    return {"num_insert_sites": num_insert_sites, "num_ta_sites": num_ta_sites, "num_ta_insert_sites": num_ta_insert_sites}
    
def get_subgraphs(graph_dir, graph_type):
    """ DEPRECATED
    Get subgraphs from the specified directory and type.

    Args:
        graph_dir (Path): Path to the directory containing graph data
        graph_type (str): Type of graph (e.g., 'case', 'control')

    Returns:
        dict: Dictionary mapping chromosome names to a list of subgraphs
    """
    subgraph_dict = {}
    for graph in (graph_dir / graph_type).iterdir():
        chrom = graph.name.split(".")[0]
        G = nx.read_graphml(graph)
        subgraphs_by_nodes = sorted(nx.connected_components(G), key=len, reverse=True)
        subgraphs = [ G.subgraph(x) for x in subgraphs_by_nodes ]

        subgraph_dict[chrom] = subgraphs
    return subgraph_dict

def get_subgraph_stats(subgraphs, graph_type, chrom, bed_df, ta_error):
    subgraph_df_list = []
    for i, subgraph in enumerate(subgraphs):
        sg_meta = {"type": graph_type, "chrom": chrom, "subgraph": i}
        sg_prop = subgraph_properties(subgraph)
        sg_ta = subgraph_TA_sites(subgraph, bed_df, ta_error)
        sg_df = pd.DataFrame((sg_meta | sg_prop | sg_ta), index=[0])  # index doesn't matter
        subgraph_df_list.append(sg_df)
    if len(subgraph_df_list) != 0:
        return pd.concat(subgraph_df_list, ignore_index=True)
    else:
        return pd.DataFrame()
 
def pcis_overlaps(case_df, control_df):
    """
    Find overlaps between pCISs in the target DataFrame and reference DataFrame.

    Args:
        target_df (pandas.DataFrame): DataFrame of target pCISs
        reference_df (pandas.DataFrame): DataFrame of reference pCISs

    Returns:
        dict: Dictionary mapping target pCIS indexes to a list of overlapping reference pCIS indexes
    """

    # find overlapping subgraphs between cases and controls
    no_overlaps = []

    # cases overlapping with controls
    case_one_overlap = []
    case_overlaps = []
    for ca in case_df.itertuples():
        overlap_list = []
        for co in control_df.itertuples():
            if (ca.min_pos <= co.max_pos) and (ca.max_pos >= co.min_pos):
                overlap_list.append(co.subgraph)
        if len(overlap_list) == 0:
            no_overlaps.append( (ca.subgraph, np.nan) )
        elif len(overlap_list) == 1:
            case_one_overlap.append( (ca.subgraph, overlap_list[0]) )
        else:
            case_overlaps.append( (ca.subgraph, overlap_list) )
            
    # controls overlapping with cases
    control_one_overlap = []
    contol_overlaps = []
    for co in control_df.itertuples():
        overlap_list = []
        for ca in case_df.itertuples():
            if (co.min_pos <= ca.max_pos) and (co.max_pos >= ca.min_pos):
                overlap_list.append(ca.subgraph)
        if len(overlap_list) == 0:
            no_overlaps.append( (np.nan, co.subgraph) )
        elif len(overlap_list) == 1:
            control_one_overlap.append( (overlap_list[0], co.subgraph) )
        else:
            contol_overlaps.append( (overlap_list, co.subgraph) )
            
    # remove single overlaps in opposite cases/controls using the multiple overlaps
    for overlap in contol_overlaps:
        for ca in overlap[0]:
            case_one_overlap.remove( (ca, overlap[1]) )
            
    for overlap in case_overlaps:
        for co in overlap[1]:
            control_one_overlap.remove( (overlap[0], co) )

    # now case and control one overlaps should be the same
    assert sorted(case_one_overlap) == sorted(control_one_overlap)

    # join all non-overlapping, overlapping once, and overlapping multiple subgraphs into one dataframe to iterate through
    # these are the pseudo Common Insertion Sites
    all_overlaps = no_overlaps + case_one_overlap
    all_overlaps.extend(contol_overlaps)
    all_overlaps.extend(case_overlaps)
    return pd.DataFrame(all_overlaps, columns=["case", "control"], dtype="object")

def compare_pcis(overlap_df, case_chrom_subgraphs, control_chrom_subgraphs, case, control, num_cases, num_controls, chrom):
    """
    Compare pCISs between overlapping subgraphs (case vs. control).

    Args:
        target_overlaps (dict): Dictionary of overlapping pCISs
        target_subgraphs (list): List of target subgraphs
        reference_subgraphs (list): List of reference subgraphs

    Returns:
        pandas.DataFrame: DataFrame containing pCIS comparison results
    """
    
    # do case-control comparison of TAs between overlapping subgraphs (pCIS):
    # - Match TA to TA site per subgraph
    # - calculate log fold changes
    # - use binomtest for significance of each TA
    # - Then use overall statistic for independent sample test between target and all references
    # which can be used for the final determination if the pCIS is now a CIS

    IS_df_list = []
    pCIS_df_list = []
    for overlap in overlap_df.itertuples():
        # get normalized read counts (CPM) for each pCIS
        case_ind = overlap.case
        control_ind = overlap.control
        
        # TODO: simplify this with functions
        # if multiple case 
        if type(case_ind) is list:
            # get single control
            control_G = control_chrom_subgraphs[control_ind]
            control_pos = [ control_G.nodes[node]['position'] for node in control_G.nodes ]
            tmp_control = pd.DataFrame([ {"control_count": control_G.nodes[node]['CPM']} for node in control_G.nodes ], index=control_pos)
            tmp_control["control_index"] = control_ind
            
            case_samples = set()
            num_control_samples = len({ x for y in [ control_G.nodes[node]["sample_IDs"] for node in control_G.nodes ] for x in y })
            
            # get multiple cases
            tmp_case_list = []
            tmp_case_pos = []
            for case_index in case_ind:
                case_G = case_chrom_subgraphs[case_index]
                case_position = [ case_G.nodes[node]['position'] for node in case_G.nodes ]
                tmp_case_pos.extend(case_position)
                tmp_case = pd.DataFrame([ {"case_count": case_G.nodes[node]['CPM']} for node in case_G.nodes ], index=case_position)
                tmp_case["case_index"] = int(case_index)
                tmp_case_list.append(tmp_case)
                case_samples = case_samples.union({ x for y in [ case_G.nodes[node]["sample_IDs"] for node in case_G.nodes ] for x in y })
            num_case_samples = len(case_samples)
            tmp_cases = pd.concat(tmp_case_list, axis=0)
            tmp = tmp_control.join(tmp_cases, how="outer")
            
            case_pos_min = min(tmp_case_pos)
            case_pos_max = max(tmp_case_pos)
            control_pos_min = min(control_pos)
            control_pos_max = max(control_pos)

            case_IS = len(tmp_cases)
            control_IS = len(tmp_control)
            
        # if mulitple control
        elif type(control_ind) is list:
            # get single case
            case_G = case_chrom_subgraphs[case_ind]
            case_pos = [ case_G.nodes[node]['position'] for node in case_G.nodes ]
            tmp_case = pd.DataFrame([ {"case_count": case_G.nodes[node]['CPM']} for node in case_G.nodes ], index=case_pos)
            tmp_case["case_index"] = case_ind

            num_case_samples = len({ x for y in [ case_G.nodes[node]["sample_IDs"] for node in case_G.nodes ] for x in y })
            control_samples = set()
            
            # get multiple controls
            tmp_control_list = []
            tmp_control_pos = []
            for control_index in control_ind:
                control_G = control_chrom_subgraphs[control_index]
                control_pos = [ control_G.nodes[node]['position'] for node in control_G.nodes ]
                tmp_control_pos.extend(control_pos)
                tmp_control = pd.DataFrame([ {"control_count": control_G.nodes[node]['CPM']} for node in control_G.nodes ], index=control_pos)
                tmp_control["control_index"] = int(control_index)
                tmp_control_list.append(tmp_control)
                control_samples = control_samples.union({ x for y in [ control_G.nodes[node]["sample_IDs"] for node in control_G.nodes ] for x in y })
            num_control_samples = len(control_samples)
            tmp_controls = pd.concat(tmp_control_list, axis=0)
            tmp = tmp_case.join(tmp_controls, how="outer")
            
            case_pos_min = min(case_pos)
            case_pos_max = max(case_pos)
            control_pos_min = min(tmp_control_pos)
            control_pos_max = max(tmp_control_pos)

            case_IS = len(tmp_case)
            control_IS = len(tmp_controls)
            
        # if just case 
        elif control_ind is None or np.isnan(control_ind):
            case_G = case_chrom_subgraphs[case_ind]
            case_pos = [ case_G.nodes[node]['position'] for node in case_G.nodes ]
            tmp_case = pd.DataFrame([ {"case_count": case_G.nodes[node]['CPM']} for node in case_G.nodes ], index=case_pos)
            tmp_case["case_index"] = case_ind
            
            tmp = tmp_case
            tmp["control_count"] = 0.0
            tmp["control_index"] = np.nan
            
            num_case_samples = len({ x for y in [ case_G.nodes[node]["sample_IDs"] for node in case_G.nodes ] for x in y })
            num_control_samples = 0
            
            case_pos_min = min(case_pos)
            case_pos_max = max(case_pos)
            control_pos_min = None
            control_pos_max = None
            
            case_IS = len(tmp)
            control_IS = 0
            
        # if just control
        elif case_ind is None or np.isnan(case_ind):
            control_G = control_chrom_subgraphs[control_ind]
            control_pos = [ control_G.nodes[node]['position'] for node in control_G.nodes ]
            tmp_control = pd.DataFrame([ {"control_count": control_G.nodes[node]['CPM']} for node in control_G.nodes ], index=control_pos)
            tmp_control["control_index"] = control_ind
            
            tmp = tmp_control
            tmp["case_count"] = 0.0
            tmp["case_index"] = np.nan
            
            num_case_samples = 0
            num_control_samples = len({ x for y in [ control_G.nodes[node]["sample_IDs"] for node in control_G.nodes ] for x in y })
        
            case_pos_min = None
            case_pos_max = None
            control_pos_min = min(control_pos)
            control_pos_max = max(control_pos)

            case_IS = 0
            control_IS = len(tmp)
            
        # if one case and one control
        elif type(case_ind) is not list and type(control_ind) is not list:
            case_G = case_chrom_subgraphs[case_ind]
            case_pos = [ case_G.nodes[node]['position'] for node in case_G.nodes ]
            tmp_case = pd.DataFrame([ {"case_count": case_G.nodes[node]['CPM']} for node in case_G.nodes ], index=case_pos)
            tmp_case["case_index"] = case_ind
            
            control_G = control_chrom_subgraphs[control_ind]
            control_pos = [ control_G.nodes[node]['position'] for node in control_G.nodes ]
            tmp_control = pd.DataFrame([ {"control_count": control_G.nodes[node]['CPM']} for node in control_G.nodes ], index=control_pos)
            tmp_control["control_index"] = control_ind
            
            tmp = tmp_case.join(tmp_control, how="outer")
            
            num_case_samples = len({ x for y in [ case_G.nodes[node]["sample_IDs"] for node in case_G.nodes ] for x in y })
            num_control_samples = len({ x for y in [ control_G.nodes[node]["sample_IDs"] for node in control_G.nodes ] for x in y })
            
            case_pos_min = min(case_pos)
            case_pos_max = max(case_pos)
            control_pos_min = min(control_pos)
            control_pos_max = max(control_pos)
            
            case_IS = len(tmp_case)
            control_IS = len(tmp_control)
            
        else:
            print(overlap)
            print("this shouldn't happen")
        
        # only fillna in case_count and control_count
        tmp["case_count"] = tmp["case_count"].fillna(0.0)
        tmp["control_count"] = tmp["control_count"].fillna(0.0)
        tmp = tmp.reset_index(drop=False).rename(columns={"index": "pos"})


        # run stats for each pCIS
        # NOTE: binomtest takes only integeres, so I'm converting the normalized read counts to the closest integers
        # get stats per TA site (only count is used)
        # used pseudo count of 1 for log fold change, and so I wanted to show the difference in binomial test and significance with this
        tmp["target_binom_pval"] = tmp.apply(lambda x: binomtest( int(x["case_count"]) + 1, int(x["case_count"] + x["control_count"]) + 1 ).pvalue, axis=1)
        tmp["target_binom_sig"] = tmp["target_binom_pval"] < 0.05
        tmp["LFC"] = tmp.apply(lambda x: np.log2((x["case_count"] + 1) / (x["control_count"] + 1)), axis=1)

        rs = ranksums(tmp["case_count"], tmp["control_count"]).pvalue
        binom = binomtest(int(tmp["case_count"].sum()) + 1, int(tmp["case_count"].sum() + tmp["control_count"].sum()) + 1, 0.5).pvalue
        # contingency table = [[a, b], [c, d]]
        #            in pCIS   not in pCIS
        # target        a           b
        # reference     c           d
        a = num_case_samples
        b = num_cases - num_case_samples
        c = num_control_samples
        d = num_controls - num_control_samples
        if a < 0 or b < 0 or c < 0 or d < 0:
            print(chrom)
            print(a, b, num_cases)
            print(c, d, num_controls)
        fi = fisher_exact([[a, b], [c, d]]).pvalue
        
        lfc = np.log2( (tmp["case_count"].sum()+1) / (tmp["control_count"].sum()+1) )
        total_IS = len(tmp)
        tmp2 = {
            "case_index": case_ind,
            "case_pos_min": case_pos_min,
            "case_pos_max": case_pos_max,
            
            "control_index": control_ind,
            "control_pos_min": control_pos_min,
            "control_pos_max": control_pos_max,
            
            "LFC": lfc,
            "ranksums": rs,
            "binomial": binom,
            "fishers_exact": fi,
        
            "total_num_samples": num_case_samples + num_control_samples,
            "case_num_samples": num_case_samples,
            "control_num_samples": num_control_samples,
            
            "total_IS": total_IS,
            "case_IS": case_IS,
            "control_IS": control_IS,
            
            "case_total_read_count": tmp["case_count"].sum(),
            "control_total_read_count": tmp["control_count"].sum(),
            }
        
        IS_df_list.append(tmp)
        pCIS_df_list.append(tmp2)
        
    IS_df = pd.concat(IS_df_list, ignore_index=True)
    IS_df["case"] = case
    IS_df["control"] = control
    IS_df["chrom"] = chrom

    pCIS_df = pd.DataFrame(pCIS_df_list)
    pCIS_df["case"] = case
    pCIS_df["control"] = control
    pCIS_df["chrom"] = chrom
    return IS_df, pCIS_df

def chrom_analysis(iter_args):
    chrom, annot_chrom_df, chrom_bed_file, args = iter_args
    graph_dir = args["graph_dir"]
    case = args["case"]
    control = args["control"]
    ta_error = args["ta_error"]
    pval_threshold = args["pval_threshold"]
    verbose = args["verbose"]
    gene_expander = args["gene_expander"]
    
    
    bed_chrom_df = pd.read_csv(chrom_bed_file, sep="\t", header=None)
    
    
    with open(graph_dir / case / chrom / "subgraphs.pickle", 'rb') as f:
        case_chrom_subgraphs = pickle.load(f)
    case_chrom_df = get_subgraph_stats(case_chrom_subgraphs, case, chrom, bed_chrom_df, ta_error)
    
    with open(graph_dir / control / chrom / "subgraphs.pickle", 'rb') as f:
        control_chrom_subgraphs = pickle.load(f)
    control_chrom_df = get_subgraph_stats(control_chrom_subgraphs, control, chrom, bed_chrom_df, ta_error)
    
    # get total samples for case and controls
    # double list comprehension https://stackoverflow.com/questions/17657720/python-list-comprehension-double-for
    case_samples = { x for y in case_chrom_df["sample_IDs"] for x in y } if len(case_chrom_df) else set()
    control_samples = { x for y in control_chrom_df["sample_IDs"] for x in y } if len(control_chrom_df) else set()
    num_cases = len(case_samples)
    num_controls = len(control_samples)
    
    # find the overlapping pCIS between case and control subgraphs
    overlap_df = pcis_overlaps(case_chrom_df, control_chrom_df)

    # get statistics for each insertion site and pCIS
    IS_df, pCIS_df = compare_pcis(overlap_df, case_chrom_subgraphs, control_chrom_subgraphs, case, control, num_cases, num_controls, chrom)
    
    # trim down annotation dataframe to just genes
    annot_chrom_genes = annot_chrom_df[annot_chrom_df["Marker Type"] == "Gene"]
    gene_names = annot_chrom_genes["Marker Symbol"].to_numpy()
    # find genes within a pCIS that includes the gene expander range
    pos_min = pCIS_df[["case_pos_min", "control_pos_min"]].min(axis=1).to_numpy().reshape(-1, 1)
    pos_max = pCIS_df[["case_pos_max", "control_pos_max"]].max(axis=1).to_numpy().reshape(-1, 1)
    gene_start = (annot_chrom_genes["genome coordinate start"] - gene_expander).to_numpy().reshape(1, -1)
    gene_end = (annot_chrom_genes["genome coordinate end"] + gene_expander).to_numpy().reshape(1, -1)
    tmp = (pos_min <= gene_end) & (pos_max >= gene_start)
    # add on genes to pCIS
    pCIS_df["genes"] = [ list(gene_names[tmp[i]]) for i in range(tmp.shape[0]) ]
    
    # # filter down to significant CIS
    # fet = pCIS_df["fishers_exact"] <= pval_threshold
    # rst = pCIS_df["ranksums"] <= pval_threshold
    # bit = pCIS_df["binomial"] <= pval_threshold
    # CIS_df = pCIS_df[ fet | rst | bit ]
    
    print(f"{chrom}\t{len(pCIS_df)}")  # /{len(pCIS_df)}")
    
    return {"is": IS_df, "pcis": pCIS_df}  # , "cis": CIS_df}

def main(args):
    """
    Main function to perform CIS analysis.
    
    Args:
        args (dict): Command-line arguments
        
    Returns:
        None
    """
    case = args["case"]
    control = args["control"]
    graph_dir = args["graph_dir"]
    
    output_res = args["output"] / f"{case}-{control}"
    output_res.mkdir(exist_ok=True)
    
    annot_df = pd.read_csv(args["gene_annot"], sep="\t")
    annot_df = annot_df[pd.notna(annot_df["genome coordinate start"])].drop("Status", axis=1)
    annot_df["chrom"] = annot_df["Chr"].apply(lambda x: f"chr{x}")
    annot_df = annot_df.sort_values("chrom")
    # TODO: what about the strand in annot_df?

    # bed_files = { file.name.split(".")[0]: pd.read_csv(file, sep="\t", header=None) for file in args["ta_dir"].iterdir() }
    bed_files = {file.name.split(".")[0]: file for file in args["ta_dir"].iterdir()}
   
    # chroms = case_df["chrom"].sort_values().unique()
    chroms = sorted([ chrom.name for chrom in (graph_dir / case).iterdir() ])
    
    # don't allow more jobs than there are chromosomes
    jobs = args["njobs"]
    num_chr = len(chroms)
    if num_chr < jobs:
        print(f"Reducing number of jobs from {jobs} to {num_chr}, since there are only {num_chr} chromosomes present.")
        jobs = len(chroms)
            
    # iter_args = tqdm([ (chrom, annot_df[annot_df["chrom"] == chrom], bed_files[chrom], args) for chrom in chroms ])
    iter_args = [ (chrom, annot_df[annot_df["chrom"] == chrom], bed_files[chrom], args) for chrom in chroms ]
    print("chrom\t# pCIS")
    with Pool(args["njobs"]) as p:
        res_dict_list = [ x for x in p.imap_unordered(chrom_analysis, iter_args) ]


    # join chromosomes results together  
    IS_list = []
    pCIS_list = []
    # CIS_list = []
    for res_dict in res_dict_list:
        IS_list.append(res_dict["is"])
        pCIS_list.append(res_dict["pcis"])
        # CIS_list.append(res_dict["cis"])

    IS_df = pd.concat(IS_list, ignore_index=True)
    pCIS_df = pd.concat(pCIS_list, ignore_index=True)
    # CIS_df = pd.concat(CIS_list, ignore_index=True)

    # save results
    IS_df.to_csv(output_res / "IS.tsv", sep="\t", index=False)
    pCIS_df.to_csv(output_res / "pCIS.tsv", sep="\t", index=False)
    # CIS_df.to_csv(output_res / "CIS.tsv", sep="\t", index=False)

if __name__ == "__main__": 
    main(load_args())
