import pickle
from pathlib import Path
from multiprocessing import Pool

from docopt import docopt
import pandas as pd
import numpy as np
import seaborn.objects as so
from seaborn import axes_style
import networkx as nx
from scipy.stats import chisquare, binomtest, ranksums, mannwhitneyu, skewtest, kurtosistest
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
    TODO: add in ability to take in multiple cases/controls?
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
     -p, --pval_threshold              p-value to exclude pCIS for significance [default: 0.05]
     -j, --jobs=N                      number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "ta_error", "jobs"]
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
    num_inserts = sum([ G.nodes[node]['counts'] for node in G.nodes ])
    subgraphs_by_nodes = sorted(nx.connected_components(G), key=len, reverse=True)
    num_subgraphs = len(subgraphs_by_nodes)
    if verbose:
        print(f"number of nodes: {nodes}")
        print(f"number of edges: {edges}")
        print(f"number of insertions: {num_inserts}")
        print(f"number of subgraphs (pCIS) {num_subgraphs}")
    return {"nodes": nodes, "edges": edges, "num_inserts": num_inserts, "num_subgraphs": num_subgraphs}
    
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
    num_inserts = sum([ G.nodes[node]['counts'] for node in G.nodes ])
    tmp_pos = sorted([ G.nodes[node]["position"] for node in G.nodes ])
    min_pos = min(tmp_pos)
    max_pos = max(tmp_pos)
    range_pos = max_pos - min_pos
    if verbose:
        print(f"number of nodes: {nodes}")
        print(f"number of edges: {edges}")
        print(f"number of insertions: {num_inserts}")
        print(f"min position: {min_pos}")
        print(f"max position: {max_pos}")
        print(f"range: {range_pos}")
    # TODO: FIXME: total nodes and range are not inclusive. see chr7
    # 7/10/23 - what does this mean...?
    return {"nodes": nodes, "edges": edges, "num_inserts": num_inserts, "min_pos": min_pos, "max_pos": max_pos, "range": range_pos}
    
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
 
def pcis_overlaps(target_df, reference_df):
    """
    Find overlaps between pCISs in the target DataFrame and reference DataFrame.

    Args:
        target_df (pandas.DataFrame): DataFrame of target pCISs
        reference_df (pandas.DataFrame): DataFrame of reference pCISs

    Returns:
        dict: Dictionary mapping target pCIS indexes to a list of overlapping reference pCIS indexes
    """
    # for all pCISs in target_df, find any overlap with reference_df (could be more than one or none)
    
    # start with each target subgraph, look through all reference subgraphs and record any overlap in the positions
    # key: index of target subgraph, value: list of reference subgraph indexes that overlap within the range of target subgraph
    overlap_dict = {}
    for tar_sg in target_df.itertuples():
        overlap_list = []
        for ref_sg in reference_df.itertuples():
            if (tar_sg.min_pos <= ref_sg.max_pos) and (tar_sg.max_pos >= ref_sg.min_pos):
                overlap_list.append(ref_sg.subgraph)
        overlap_dict[tar_sg.subgraph] = overlap_list
    # TODO: would it be easier/more efficient to compute this by whole dataframe?
    # this is the fully coded out logic, but now can I make it better and improve my pandas skills?
    # I can use chrom_overlaps as the ground truth as well.
    return overlap_dict

def compare_pcis(target_overlaps, target_subgraphs, reference_subgraphs, target, reference, chrom):
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
    TA_df_list = []
    overall_df_list = []

    for tar_ind, ref_inds in target_overlaps.items():
        tar_G = target_subgraphs[tar_ind]
        tar_pos = [ tar_G.nodes[node]['position'] for node in tar_G.nodes ]
        tmp_tar = pd.DataFrame([ {"target_count": tar_G.nodes[node]['counts']} for node in tar_G.nodes ], index=tar_pos)        
        if len(ref_inds) == 0:
            tmp = tmp_tar
            tmp["reference_count"] = 0
            tmp["reference_index"] = np.nan
        else:
            tmp_control_list = []
            for ref_ind in ref_inds:
                ref_G = reference_subgraphs[ref_ind]
                ref_pos = [ ref_G.nodes[node]['position'] for node in ref_G.nodes ]
                tmp_control = pd.DataFrame([ {"reference_count": ref_G.nodes[node]['counts']} for node in ref_G.nodes ], index=ref_pos)
                tmp_control["reference_index"] = int(ref_ind)
                tmp_control_list.append(tmp_control)
            tmp_control = pd.concat(tmp_control_list, axis=0)
            tmp = tmp_tar.join(tmp_control, how="outer")
            tmp = tmp.fillna(0).astype(int)
            
        tmp = tmp.reset_index(drop=False).rename(columns={"index": "pos"})
        tmp["target_index"] = tar_ind
        
        # get stats per TA site (only count is used)
        tmp["target_binom_pval"] = tmp.apply(lambda x: binomtest(int(x["target_count"]), int(x["target_count"]) + int(x["reference_count"])).pvalue, axis=1)
        tmp["target_binom_sig"] = tmp["target_binom_pval"] < 0.05
        tmp["LFC"] = tmp.apply(lambda x: np.log2((x["target_count"]+1) / (x["reference_count"]+1)), axis=1)
        # used pseudo count of 1 for log fold change, and so I wanted to show the difference in binomial test and significance with this
        tmp["p_target_binom_pval"] = tmp.apply(lambda x: binomtest(x["target_count"]+1, (x["target_count"]+1) + (x["reference_count"]+1)).pvalue, axis=1)
        tmp["p_target_binom_sig"] = tmp["p_target_binom_pval"] < 0.05
        
        # overall test stat for independence. use genomic positions. Total samples are each position times counts.
        # ex.) pos: 1001 and count: 3 is [1001, 1001, 1001]
        target_overall = []
        reference_overall = []
        for row in tmp.itertuples():
            for pos_tmp in [row.pos] * int(row.target_count):
                target_overall.append(int(pos_tmp))
            for pos_tmp in [row.pos] * int(row.reference_count):
                reference_overall.append(int(pos_tmp))

        mwu = mannwhitneyu(target_overall, reference_overall).pvalue if len(reference_overall) != 0 else np.nan
        rs = ranksums(target_overall, reference_overall).pvalue if len(reference_overall) != 0 else np.nan
        # case_skewtest = skewtest(target_overall).pvalue if len(target_overall) >= 8 else np.nan
        # case_kurtosistest = kurtosistest(target_overall).pvalue if len(target_overall) >= 20 else np.nan
        # control_skewtest = skewtest(reference_overall).pvalue if len(reference_overall) >= 8 else np.nan
        # control_kurtosistest = kurtosistest(reference_overall).pvalue if len(reference_overall) >= 20 else np.nan
        total_TA = len(tmp)
        TA_sig = tmp["target_binom_sig"].sum()
                            
        # Other stats: Kurtosis, skewness, etc. # "stat type": ["statistic", "pvalue"],
        tmp2 = pd.DataFrame({
            "target_index": [tar_ind],
            "reference_index": [tmp["reference_index"].values[0]],
            "target_pos_min": [min(tar_pos)],
            "target_pos_max": [max(tar_pos)],
            "reference_pos_min": [int(min(ref_pos))] if len(ref_inds) != 0 else [None],
            "reference_pos_max": [int(max(ref_pos))] if len(ref_inds) != 0 else [None],
            "mannwhitneyu": [mwu],
            "ranksums": [rs], 
            #  "case-skewtest": [case_skewtest],
            #  "case-kurtosistest": [case_kurtosistest],
            #  "control-skewtest": [control_skewtest],
            #  "control-kurtosistest": [control_kurtosistest],
            "total_TA": [total_TA],
            "TA_sig": [TA_sig],
            })
        TA_df_list.append(tmp)
        overall_df_list.append(tmp2)
    
    TA_df = pd.concat(TA_df_list, ignore_index=True)
    TA_df["target"] = target
    TA_df["reference"] = reference
    TA_df["chrom"] = chrom
    overall_df = pd.concat(overall_df_list, ignore_index=True)
    overall_df["sig_ratio"] = overall_df["TA_sig"] / overall_df["total_TA"]
    overall_df["target"] = target
    overall_df["reference"] = reference
    overall_df["chrom"] = chrom
    return TA_df, overall_df

def pcis_to_cis(overall_df, threshold):
    """
    Convert pCISs to CISs based on statistical significance.

    Args:
        overall_df (pandas.DataFrame): DataFrame containing overall pCIS statistics
        threshold (float): Significance threshold for p-values

    Returns:
        pandas.DataFrame: DataFrame containing significant CISs
    """
    # find pcis with significant pvalue that is less than the given threshold
    sig_df = overall_df[ (overall_df["mannwhitneyu"] < threshold) & (overall_df["ranksums"] < threshold) ]
    # test stat below threshold OR ratio that is not 0 and there are more than 1 sig tA
    nan_df =  overall_df[ pd.isna(overall_df["mannwhitneyu"]) & pd.isna(overall_df["ranksums"]) ]
    nan_sig_df = nan_df[ (nan_df["sig_ratio"] != 0) & (nan_df["TA_sig"] > 1) ]
    all_sig_df = pd.concat([sig_df, nan_sig_df]).reset_index(drop=True)
    return all_sig_df

def cis_annotate(target_sig_df, annotated_df, gene_expander=50000):
    """
    Annotate CISs with gene markers.

    Args:
        target_sig_df (pandas.DataFrame): DataFrame containing significant CISs
        annotated_df (pandas.DataFrame): DataFrame containing gene annotations

    Returns:
        pandas.DataFrame: DataFrame containing annotated CISs
    """
    # get genes markers to each CIS
    gene_list = []
    for row in target_sig_df.itertuples():
        sub_gene_list = []
        for gene in annotated_df.itertuples():
            try:
                sub3 = row.reference_pos_min <= (gene._5 + gene_expander)
                sub4 = row.reference_pos_max >= (gene._4 - gene_expander)
            except:
                sub3 = None
                sub4 = None
                
            if (row.target_pos_min <= (gene._5 + gene_expander)) and (row.target_pos_max >= (gene._4 - gene_expander)):
                sub_gene_list.append(pd.DataFrame({
                    "type": ["target"],
                    "type_name": [row.target],
                    "type_index": [int(row.target_index)],
                    "chrom": [row.chrom],
                    "marker_symbol": [gene._7],
                    "marker_name": [gene._8],
                    "marker_type": [gene._9],
                    "marker_feature_type": [gene._10],
                    "marker_annot_index": [gene.Index],
                    "genome coordinate start": [gene._4],
                    "genome coordinate end": [gene._5],
                    "genome coordinate expander": [gene_expander]
                    }))
            elif sub3 and sub4:
                sub_gene_list.append(pd.DataFrame({
                    "type": ["reference"],
                    "type_name": [row.reference],
                    "type_index": [int(row.reference_index)],
                    "chrom": [row.chrom],
                    "marker_symbol": [gene._7],
                    "marker_name": [gene._8],
                    "marker_type": [gene._9],
                    "marker_feature_type": [gene._10],
                    "marker_annot_index": [gene.Index],
                    "genome coordinate start": [gene._4],
                    "genome coordinate end": [gene._5],
                    "genome coordinate expander": [gene_expander]
                    }))
                
        # it is possible that there were no annotations found
        if len(sub_gene_list) == 0:
            sub_gene_list.append(pd.DataFrame({
                "type": ["target"],
                "type_name": [row.target],
                "type_index": [int(row.target_index)],
                "chrom": [row.chrom],
                "marker_symbol": [None],
                "marker_name": [None],
                "marker_type": [None],
                "marker_feature_type": [None],
                "marker_annot_index": [None],
                "genome coordinate start": [None],
                "genome coordinate end": [None],
                "genome coordinate expander": [gene_expander]
                }))
            
        gene_list.extend(sub_gene_list)
    return pd.concat(gene_list, ignore_index=True)

def chrom_analysis(iter_args):
    chrom, annot_chrom_df, chrom_bed_file, args = iter_args
    graph_dir = args["graph_dir"]
    case = args["case"]
    control = args["control"]
    ta_error = args["ta_error"]
    pval_threshold = args["pval_threshold"]
    verbose = args["verbose"]
    gene_expander = 50000  # TODO: add to input args
    
    
    bed_chrom_df = pd.read_csv(chrom_bed_file, sep="\t", header=None)
    
    with open(graph_dir / case / chrom / "subgraphs.pickle", 'rb') as f:
        case_chrom_subgraphs = pickle.load(f)
    case_chrom_df = get_subgraph_stats(case_chrom_subgraphs, case, chrom, bed_chrom_df, ta_error)
    
    with open(graph_dir / control / chrom / "subgraphs.pickle", 'rb') as f:
        control_chrom_subgraphs = pickle.load(f)
    control_chrom_df = get_subgraph_stats(control_chrom_subgraphs, control, chrom, bed_chrom_df, ta_error)

    
    # cases as the target
    case_overlaps = pcis_overlaps(case_chrom_df, control_chrom_df)
    if not case_overlaps:  # if empty
        case_features, case_TA_df, case_overall_df, case_sig_df = None, None, None, None
    else:
        case_TA_df, case_overall_df = compare_pcis(case_overlaps, case_chrom_subgraphs, control_chrom_subgraphs, case, control, chrom)
        case_sig_df = pcis_to_cis(case_overall_df, pval_threshold)
        if len(case_sig_df) != 0:
            case_features = cis_annotate(case_sig_df, annot_chrom_df, gene_expander)
        else:
            case_features = None
    
    # controls as the target
    control_overlaps = pcis_overlaps(control_chrom_df, case_chrom_df)
    if not control_overlaps:  # if empty
        control_features, control_TA_df, control_overall_df, control_sig_df = None, None, None, None
    else:
        control_TA_df, control_overall_df = compare_pcis(control_overlaps, control_chrom_subgraphs, case_chrom_subgraphs, control, case, chrom)
        control_sig_df = pcis_to_cis(control_overall_df, pval_threshold)
        if len(control_sig_df) != 0:
            control_features = cis_annotate(control_sig_df, annot_chrom_df, gene_expander)
        else:
            control_features = None
    
    if case_features is not None or control_features is not None:
        genomic_features_df = pd.concat([case_features, control_features], ignore_index=True)
        if verbose:
            print(f"""{chrom}\tsig. genomic features: {genomic_features_df["marker_symbol"].unique().shape[0]}/{annot_chrom_df["Marker Symbol"].unique().shape[0]}""")
    else:
        genomic_features_df = None
        if verbose:
            print(f"{chrom}\tno sig. genomic features found")

    ta_df = pd.concat([case_TA_df, control_TA_df], ignore_index=True)
    overall_df = pd.concat([case_overall_df, control_overall_df], ignore_index=True)
    sig_df = pd.concat([case_sig_df, control_sig_df], ignore_index=True)
    graph_chrom_df = pd.concat([case_chrom_df, control_chrom_df], ignore_index=True)
    
    return {"ta": ta_df, "overall": overall_df, "sig": sig_df, "genomic_features": genomic_features_df, "graph_stats": graph_chrom_df}

def volcano_plot(data, lfc, pval, threshold=0.05):
    """
    Create a volcano plot to visualize p-value and log fold change.

    Args:
        data (pandas.DataFrame): Data for the plot
        lfc (str): Column name for log fold change values
        pval (str): Column name for p-values
        threshold (float, optional): Significance threshold for p-values. Default is 0.05.

    Returns:
        seaborn.axisgrid.FacetGrid: Volcano plot visualization
    """
    thres = np.log10(threshold) * -1
    data[pval] = np.log10(data[pval]) * -1
    g = (
        so.Plot(data, x=lfc, y=pval, pointsize=pval)
        # .add(so.Dots(), so.Jitter(1))
        .add(so.Dots(color="grey"), so.Jitter(1), data=data.query(f"{pval} < {thres}"))
        .add(so.Dots(color="blue"), so.Jitter(1), data=data.query(f"{pval} >= {thres}"))
        .scale(y="log")
    )
    return g
                
def aak1(gene_list):
    """
    Analyze the gene list for presence of specific genes and summary statistics.

    Args:
        gene_list (pandas.DataFrame): DataFrame containing gene annotations

    Returns:
        None
    """
    print(f'Number CIS without annotation: {pd.isna(gene_list["annot_index"]).sum()}')
    a_genes = gene_list["marker_symbol"].unique()
    t_genes = gene_list[gene_list["type"] == "target"]["marker_symbol"].unique()
    r_genes = gene_list[gene_list["type"] == "reference"]["marker_symbol"].unique()
    print(f"all genes: {len(a_genes)}")
    print(f"target genes: {len(t_genes)}")
    print(f"reference genes: {len(r_genes)}")
    print(f'Aak1? {"Aak1" in gene_list["marker_symbol"].unique()}')
            

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
    
    output = args["output"] / f"{case}-{control}"
    output.mkdir(exist_ok=True)
    
    annot_df = pd.read_csv(args["gene_annot"], sep="\t")
    annot_df = annot_df[pd.notna(annot_df["genome coordinate start"])].drop("Status", axis=1)
    annot_df["chrom"] = annot_df["Chr"].apply(lambda x: f"chr{x}")
    annot_df = annot_df.sort_values("chrom")
    # TODO: what about the strand in annot_df?

    # bed_files = { file.name.split(".")[0]: pd.read_csv(file, sep="\t", header=None) for file in args["ta_dir"].iterdir() }
    bed_files = {file.name.split(".")[0]: file for file in args["ta_dir"].iterdir()}
   
    # chroms = case_df["chrom"].sort_values().unique()
    chroms = sorted([ chrom.name for chrom in (args["graph_dir"] / case).iterdir() ])
    
    
    # iter_args = tqdm([ (chrom, annot_df[annot_df["chrom"] == chrom], bed_files[chrom], args) for chrom in chroms ])
    iter_args = [ (chrom, annot_df[annot_df["chrom"] == chrom], bed_files[chrom], args) for chrom in chroms ]
    with Pool(args["jobs"]) as p:
        res_dict_list = [ x for x in p.imap_unordered(chrom_analysis, iter_args) ]
        
    # save data  
    ta_list = []
    overall_list = []
    sig_list = []
    genomic_features_list = []
    graphs_stats = []
    for res_dict in res_dict_list:
        ta_list.append(res_dict["ta"])
        overall_list.append(res_dict["overall"])
        sig_list.append(res_dict["sig"])
        genomic_features_list.append(res_dict["genomic_features"])
        graphs_stats.append(res_dict["graph_stats"])


    TA_df = pd.concat(ta_list, ignore_index=True)
    TA_df.to_csv(output / "TA.tsv", sep="\t", index=False)

    overall_df = pd.concat(overall_list, ignore_index=True)
    overall_df.to_csv(output / "overall.tsv", sep="\t", index=False)

    sig_df = pd.concat(sig_list, ignore_index=True)
    sig_df.to_csv(output / "sig.tsv", sep="\t", index=False)

    genomic_features_df = pd.concat(genomic_features_list, ignore_index=True)
    genomic_features_df.to_csv(output / "genomic_features.tsv", sep="\t", index=False)

    graph_stats_df = pd.concat(graphs_stats, ignore_index=True)
    graph_stats_df.to_csv(output / "graph_stats.tsv", sep="\t", index=False)
    
if __name__ == "__main__": 
    main(load_args())
