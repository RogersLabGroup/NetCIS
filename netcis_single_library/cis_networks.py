import pickle
from pathlib import Path
from multiprocessing import Pool

from docopt import docopt
import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import binomtest, ranksums, fisher_exact


def load_args() -> dict:
    doc = """  
    Generate common insertion sites (CIS) by finding overlapping pseudo-CIS (pCIS) between cases and controls.
    CIS can still be generated even if there is no overlap found for a pCIS or if it overlaps with multiple pCIS.

    Usage: 
        cis_networks.py --output_prefix DIR --case STR --control STR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-CIS" appended to it
     -a, --case=STR                    treatment type value to use as case
     -b, --control=STR                 treatment type value to use as control
     
    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information if available using 0, 1 or 2 [default: 0]
     -j, --njobs=N                     number of processes to run [default: 1]
    """
    
    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "njobs"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])
        
    new_args["graph_dir"] = Path(new_args["output_prefix"] + "-graphs")
    new_args["output"] = Path(new_args["output_prefix"] + "-CIS")
    new_args["output"].mkdir(exist_ok=True)
    
    if new_args["verbose"] > 1:
        print(new_args)

    return new_args
    
def subgraph_properties(G):
    # Calculate properties of a subgraph.

    nodes = G.number_of_nodes()
    edges = G.number_of_edges()
    norm_num_inserts = sum([ G.nodes[node]['CPM'] for node in G.nodes ])
    tmp_pos = sorted([ G.nodes[node]["position"] for node in G.nodes ])
    min_pos = min(tmp_pos)
    max_pos = max(tmp_pos)
    range_pos = max_pos - min_pos
    
    sample_IDs = { x for y in [ G.nodes[node]["sample_IDs"] for node in G.nodes ] for x in y }
    num_unique_samples = len(sample_IDs)
        
    return {"nodes": nodes, 
            "edges": edges, 
            "norm_num_inserts": norm_num_inserts, 
            "min_pos": min_pos, 
            "max_pos": max_pos, 
            "range": range_pos,
            "sample_IDs": [sample_IDs],
            "num_unique_samples": num_unique_samples,
            }
    
def subgraph_TA_sites(G, bed, ta_error):
    # DEPRECATED but could be useful in the future if someone cares about finding TA sites per subgraph...?
    # Calculate properties related to TA sites in a subgraph.
    num_insert_sites = G.number_of_nodes()
    
    tmp_pos = sorted([ G.nodes[node]["position"] for node in G.nodes ])
    ta_sites = bed[(bed[1] > min(tmp_pos)) & (bed[2] < max(tmp_pos))]
    num_ta_sites = len(ta_sites)
    
    arr1 = np.array(tmp_pos).reshape(-1, 1)
    arr2 = ta_sites[1].to_numpy().reshape(-1, 1)
    arr3 = ta_sites[1].to_numpy().reshape(-1, 1)
    ta_inserts = (arr1 >= (arr2.T - ta_error)) & (arr1 <= (arr3.T + ta_error))
    num_ta_insert_sites = ta_inserts.any(axis=1).sum()
    
    return {"num_insert_sites": num_insert_sites, "num_ta_sites": num_ta_sites, "num_ta_insert_sites": num_ta_insert_sites}

def subgraph_stats(subgraphs, graph_type, chrom):
    subgraph_df_list = []
    for i, subgraph in enumerate(subgraphs):
        sg_meta = {"type": graph_type, "chrom": chrom, "subgraph": i}
        sg_prop = subgraph_properties(subgraph)
        # sg_ta = subgraph_TA_sites(subgraph, bed_df, ta_error)  # NOTE: 2/29/24 - removing bed files and TA sites as part of the stats per subgraph
        sg_df = pd.DataFrame((sg_meta | sg_prop), index=[0])  # index doesn't matter, it just needs one
        subgraph_df_list.append(sg_df)
    if len(subgraph_df_list) != 0:
        return pd.concat(subgraph_df_list, ignore_index=True)
    else:
        return pd.DataFrame()
 
def pcis_overlaps(case_df, control_df):
    # Find overlaps between pCISs in the target DataFrame and reference DataFrame.
    
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

def pcis_to_cis_with_stats(overlap_df, case_chrom_subgraphs, control_chrom_subgraphs, case, control, num_cases, num_controls, chrom):
    """Compare pCISs between overlapping subgraphs (case vs. control).
        
        do case-control comparison of TAs between overlapping subgraphs (pCIS):
        - Match TA to TA site per subgraph
        - calculate log fold changes
        - use binomtest for significance of each TA
        - Then use overall statistic for independent sample test between target and all references
        which can be used for the final determination if the pCIS is now a CIS
    """
    
    # create the final CIS by grouping the insertion sites (IS) from each pCIS overlap based on the conditions of the overlap
    IS_df_list = []
    CIS_df_list = []
    for overlap in overlap_df.itertuples():
        case_ind = overlap.case
        control_ind = overlap.control
        
        # TODO: refactor this with functions to be neater and shorter
        
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
            tmp_IS = tmp_control.join(tmp_cases, how="outer")
            
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
            tmp_IS = tmp_case.join(tmp_controls, how="outer")
            
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
            
            tmp_IS = tmp_case
            tmp_IS["control_count"] = 0.0
            tmp_IS["control_index"] = np.nan
            
            num_case_samples = len({ x for y in [ case_G.nodes[node]["sample_IDs"] for node in case_G.nodes ] for x in y })
            num_control_samples = 0
            
            case_pos_min = min(case_pos)
            case_pos_max = max(case_pos)
            control_pos_min = None
            control_pos_max = None
            
            case_IS = len(tmp_IS)
            control_IS = 0
            
        # if just control
        elif case_ind is None or np.isnan(case_ind):
            control_G = control_chrom_subgraphs[control_ind]
            control_pos = [ control_G.nodes[node]['position'] for node in control_G.nodes ]
            tmp_control = pd.DataFrame([ {"control_count": control_G.nodes[node]['CPM']} for node in control_G.nodes ], index=control_pos)
            tmp_control["control_index"] = control_ind
            
            tmp_IS = tmp_control
            tmp_IS["case_count"] = 0.0
            tmp_IS["case_index"] = np.nan
            
            num_case_samples = 0
            num_control_samples = len({ x for y in [ control_G.nodes[node]["sample_IDs"] for node in control_G.nodes ] for x in y })
        
            case_pos_min = None
            case_pos_max = None
            control_pos_min = min(control_pos)
            control_pos_max = max(control_pos)

            case_IS = 0
            control_IS = len(tmp_IS)
            
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
            
            tmp_IS = tmp_case.join(tmp_control, how="outer")
            
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
            assert False
        
        # if an insertion site did not appear in cases or controls, fill it in with 0
        tmp_IS["case_count"] = tmp_IS["case_count"].fillna(0.0)
        tmp_IS["control_count"] = tmp_IS["control_count"].fillna(0.0)
        tmp_IS = tmp_IS.reset_index(drop=False).rename(columns={"index": "pos"})

        # get stats for each CIS
        # binomtest takes only integeres, so I'm converting the normalized read counts to the closest integers
        # and pseudo counts of 1 are added to log fold change and binomial tests as having a 0 for either could cause errors
        lfc = np.log2( (tmp_IS["case_count"].sum() + 1) / (tmp_IS["control_count"].sum() + 1) )
        binom = binomtest(int(tmp_IS["case_count"].sum()) + 1, int(tmp_IS["case_count"].sum() + tmp_IS["control_count"].sum()) + 1, 0.5).pvalue
        rs = ranksums(tmp_IS["case_count"], tmp_IS["control_count"]).pvalue
        # contingency table = [[a, b], [c, d]]
        #                 in CIS     not in CIS
        # total cases        a           b
        # total controls     c           d
        a = num_case_samples
        b = num_cases - num_case_samples
        c = num_control_samples
        d = num_controls - num_control_samples
        table = [[a, b], [c, d]]
        fi = fisher_exact(table).pvalue
        tmp_CIS = {
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
            
            "total_IS": len(tmp_IS),
            "case_IS": case_IS,
            "control_IS": control_IS,
            
            "total_read_count": tmp_IS["case_count"].sum() + tmp_IS["control_count"].sum(),
            "case_total_read_count": tmp_IS["case_count"].sum(),
            "control_total_read_count": tmp_IS["control_count"].sum(),
            }
        
        IS_df_list.append(tmp_IS)
        CIS_df_list.append(tmp_CIS)
        
    IS_df = pd.concat(IS_df_list, ignore_index=True)
    IS_df["case"] = case
    IS_df["control"] = control
    IS_df["chrom"] = chrom

    CIS_df = pd.DataFrame(CIS_df_list)
    CIS_df["case"] = case
    CIS_df["control"] = control
    CIS_df["chrom"] = chrom
    return IS_df, CIS_df

def run_per_chrom(iter_args):
    chrom, args = iter_args
    graph_dir = args["graph_dir"]
    case = args["case"]
    control = args["control"]
    verbose = args["verbose"]
    
    with open(graph_dir / case / chrom / "subgraphs.pickle", 'rb') as f:
        case_chrom_subgraphs = pickle.load(f)
    case_chrom_df = subgraph_stats(case_chrom_subgraphs, case, chrom)
    
    with open(graph_dir / control / chrom / "subgraphs.pickle", 'rb') as f:
        control_chrom_subgraphs = pickle.load(f)
    control_chrom_df = subgraph_stats(control_chrom_subgraphs, control, chrom)
    
    # get total samples for case and controls
    # double list comprehension https://stackoverflow.com/questions/17657720/python-list-comprehension-double-for
    case_samples = { x for y in case_chrom_df["sample_IDs"] for x in y } if len(case_chrom_df) else set()
    control_samples = { x for y in control_chrom_df["sample_IDs"] for x in y } if len(control_chrom_df) else set()
    num_cases = len(case_samples)
    num_controls = len(control_samples)
    
    # find overlapping pCIS between case and control subgraphs
    overlap_df = pcis_overlaps(case_chrom_df, control_chrom_df)

    # get statistics for insertion sites and CIS. This is where pCIS -> CIS
    IS_df, CIS_df = pcis_to_cis_with_stats(overlap_df, case_chrom_subgraphs, control_chrom_subgraphs, case, control, num_cases, num_controls, chrom)
       
    if verbose:
        print(f"{chrom}\t{len(CIS_df)}")
    
    return {"is": IS_df, "cis": CIS_df}

def main(args):
    case = args["case"]
    control = args["control"]
    graph_dir = args["graph_dir"]
    
    output_res = args["output"] / f"{case}-{control}"
    output_res.mkdir(exist_ok=True)

    # bed_files = {file.name.split(".")[0]: file for file in args["ta_dir"].iterdir()}
   
    chroms = sorted([ chrom.name for chrom in (graph_dir / case).iterdir() ])
    
    # don't allow more jobs than there are chromosomes
    jobs = args["njobs"]
    num_chr = len(chroms)
    if num_chr < jobs:
        print(f"ATTENTION: Reducing number of jobs from {jobs} to {num_chr}, since there are only {num_chr} chromosomes present.")
        jobs = len(chroms)
            
    iter_args = [ (chrom, args) for chrom in chroms ]
    if args["verbose"]:
        print("chrom\t# CIS")
    with Pool(args["njobs"]) as p:
        res_dict_list = [ x for x in p.imap_unordered(run_per_chrom, iter_args) ]

    # join chromosome results together and save results 
    IS_list = [ res_dict["is"] for res_dict in res_dict_list ]
    IS_df = pd.concat(IS_list, ignore_index=True)
    IS_df.to_csv(output_res / "IS.tsv", sep="\t", index=False)
    
    CIS_list = [ res_dict["cis"] for res_dict in res_dict_list ]
    CIS_df = pd.concat(CIS_list, ignore_index=True)
    CIS_df.to_csv(output_res / "CIS.tsv", sep="\t", index=False) 
    
if __name__ == "__main__": 
    main(load_args())
