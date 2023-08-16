from pathlib import Path

from docopt import docopt
import pandas as pd
import numpy as np
import seaborn.objects as so
from seaborn import axes_style
import networkx as nx
from scipy.stats import chisquare, binomtest, ranksums, mannwhitneyu, skewtest, kurtosistest


def load_args() -> dict:
    doc = """  
    Generate common insertion sites (CIS) using a network (graph) based approach. 
    The only parameter that will change how a CIS is generated can be control with --threshold
    
    Usage: analysis.py --output_prefix DIR --ta_dir DIR --gene_annot FILE [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-analysis" appended to it
     -b, --ta_dir=DIR                  directory that contains the TA locations for each chromosome in bed format
     -g, --gene_annot=FILE             MGI's mouse menetic markers excluding withdrawn genes

    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   (TODO: how to allow --verbose meaning 1 as well as supplying value?) print more verbose information using 0, 1 or 2 [default: 0]
     -t, --ta_error=N                  how many bases to expand the search for a TA site at each insertion [default: 5]
     -j, --jobs=N                      number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "ta_error", "jobs"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])

    new_args["graph_dir"] = Path(new_args["output_prefix"] + "-graphs")
    new_args["ta_dir"] = Path(new_args["ta_dir"])
    new_args["gene_annot"] = Path(new_args["gene_annot"])
    new_args["output"] = Path(new_args["output_prefix"] + "-analysis")
    new_args["output"].mkdir(exist_ok=True)

    return new_args

def graph_properties(G, verbose=0):
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
    subgraph_dict = {}
    for graph in (graph_dir / graph_type).iterdir():
        chrom = graph.name.split(".")[0]
        G = nx.read_graphml(graph)
        subgraphs_by_nodes = sorted(nx.connected_components(G), key=len, reverse=True)
        subgraph_dict[chrom] = [ G.subgraph(x) for x in subgraphs_by_nodes ]
    return subgraph_dict

def get_subgraph_stats(subgraph_chroms, graph_type, bed_files, ta_error):
    subgraph_df_list = []
    for chrom, subgraphs in subgraph_chroms.items():
        for i, subgraph in enumerate(subgraphs):
            sg_meta = {"type": graph_type, "chrom": chrom, "subgraph": i}
            sg_prop = subgraph_properties(subgraph)
            sg_ta = subgraph_TA_sites(subgraph, bed_files[chrom], ta_error)
            sg_df = pd.DataFrame((sg_meta | sg_prop | sg_ta))
            subgraph_df_list.append(sg_df)
    return pd.concat(subgraph_df_list, ignore_index=True)
 
def pcis_overlaps(target_df, reference_df):
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

def compare_pcis(target_overlaps, target_subgraphs, reference_subgraphs):
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
        tmp_tar = pd.DataFrame([ {"target": tar_G.nodes[node]['counts']} for node in tar_G.nodes ], index=tar_pos)        
        if len(ref_inds) == 0:
            tmp = tmp_tar
            tmp["reference"] = 0
            tmp["reference_index"] = np.nan
        else:
            for ref_ind in ref_inds:
                ref_G = reference_subgraphs[ref_ind]
                ref_pos = [ ref_G.nodes[node]['position'] for node in ref_G.nodes ]
                tmp_control = pd.DataFrame([ {"reference": ref_G.nodes[node]['counts']} for node in ref_G.nodes ], index=ref_pos)                
            # get union of all insertion sites and make it into a df
            tmp = tmp_tar.join(tmp_control, how="outer")
            tmp = tmp.fillna(0).astype(int)
            tmp["reference_index"] = "-".join([str(x) for x in ref_inds])
            
        tmp = tmp.reset_index(drop=False).rename(columns={"index": "pos"})
        tmp["target_index"] = tar_ind
        
        # get stats per TA site (only count is used)
        tmp["target_binom_pval"] = tmp.apply(lambda x: binomtest(int(x["target"]), int(x["target"]) + int(x["reference"])).pvalue, axis=1)
        tmp["target_binom_sig"] = tmp["target_binom_pval"] < 0.05
        tmp["LFC"] = tmp.apply(lambda x: np.log2((x["target"]+1) / (x["reference"]+1)), axis=1)
        # used pseudo count of 1 for log fold change, and so I wanted to show the difference in binomial test and significance with this
        tmp["p_target_binom_pval"] = tmp.apply(lambda x: binomtest(x["target"]+1, (x["target"]+1) + (x["reference"]+1)).pvalue, axis=1)
        tmp["p_target_binom_sig"] = tmp["p_target_binom_pval"] < 0.05
        
        # overall test stat for independence. use genomic positions. Total samples are each position times counts.
        # ex.) pos: 1001 and count: 3 is [1001, 1001, 1001]
        target_overall = []
        reference_overall = []
        for row in tmp.itertuples():
            for pos_tmp in [row.pos] * int(row.target):
                target_overall.append(int(pos_tmp))
            for pos_tmp in [row.pos] * int(row.reference):
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
            "reference_pos_min": [min(ref_pos)] if len(ref_inds) != 0 else [None],
            "reference_pos_max": [max(ref_pos)] if len(ref_inds) != 0 else [None],
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
    overall_df = pd.concat(overall_df_list, ignore_index=True)
    overall_df["sig_ratio"] = overall_df["TA_sig"] / overall_df["total_TA"]
    return TA_df, overall_df

def pcis_to_cis(overall_df, threshold):
    # find pcis with significant pvalue that is less than the given threshold
    sig_df = overall_df[ (overall_df["mannwhitneyu"] < threshold) & (overall_df["ranksums"] < threshold) ]
    # test stat below threshold OR ratio that is not 0 and there are more than 1 sig tA
    nan_df =  overall_df[ pd.isna(overall_df["mannwhitneyu"]) & pd.isna(overall_df["ranksums"]) ]
    nan_sig_df = nan_df[ (nan_df["sig_ratio"] != 0) & (nan_df["TA_sig"] > 1) ]
    all_sig_df = pd.concat([sig_df, nan_sig_df]).reset_index(drop=True)
    return all_sig_df

def cis_annotate(target_sig_df, annotated_df):
    # get genes markers to each CIS
    gene_list = []
    for row in target_sig_df.itertuples():
        sub_gene_list = []
        for gene in annotated_df.itertuples():
            try:
                sub3 = row.reference_pos_min <= gene._5
                sub4 = row.reference_pos_max >= gene._4
            except:
                sub3 = None
                sub4 = None
                
            if (row.target_pos_min <= gene._5) and (row.target_pos_max >= gene._4):
                sub_gene_list.append(pd.DataFrame({
                    "type": ["target"],
                    "type_index": [int(row.target_index)],
                    "marker_symbol": [gene._7],
                    "marker_name": [gene._8],
                    "marker_type": [gene._9],
                    "marker_feature_type": [gene._10],
                    "marker_annot_index": [gene.Index],
                    }))
            elif sub3 and sub4:
                sub_gene_list.append(pd.DataFrame({
                    "type": ["reference"],
                    "type_index": [int(row.reference_index)],
                    "marker_symbol": [gene._7],
                    "marker_name": [gene._8],
                    "marker_type": [gene._9],
                    "marker_feature_type": [gene._10],
                    "marker_annot_index": [gene.Index],
                    }))
                
        # it is possible that there were no annotations found
        if len(sub_gene_list) == 0:
            sub_gene_list.append(pd.DataFrame({
                "type": ["target"],
                "type_index": [int(row.target_index)],
                "marker_symbol": [None],
                "marker_name": [None],
                "marker_type": [None],
                "marker_feature_type": [None],
                "marker_annot_index": [None],
                }))
            
        gene_list.extend(sub_gene_list)
    return pd.concat(gene_list, ignore_index=True)
      
def volcano_plot(data, lfc, pval, threshold=0.05):
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
    print(f'Number CIS without annotation: {pd.isna(gene_list["annot_index"]).sum()}')
    a_genes = gene_list["marker_symbol"].unique()
    t_genes = gene_list[gene_list["type"] == "target"]["marker_symbol"].unique()
    r_genes = gene_list[gene_list["type"] == "reference"]["marker_symbol"].unique()
    print(f"all genes: {len(a_genes)}")
    print(f"target genes: {len(t_genes)}")
    print(f"reference genes: {len(r_genes)}")
    print(f'Aak1? {"Aak1" in gene_list["marker_symbol"].unique()}')
            

def main(args):
    verbose = args["verbose"]
    pval_threshold = 0.05
    
    annot_df = pd.read_csv(args["gene_annot"], sep="\t")
    annot_df = annot_df[pd.notna(annot_df["genome coordinate start"])].drop("Status", axis=1)
    annot_df["chrom"] = annot_df["Chr"].apply(lambda x: f"chr{x}")
    annot_df = annot_df.sort_values("chrom")
    # TODO: what about the strand in annot_df?

    bed_files = { file.name.split(".")[0]: pd.read_csv(file, sep="\t", header=None) for file in args["ta_dir"].iterdir() }
    
    case_subgraph_dict = get_subgraphs(args["graph_dir"], "case")
    # case_df = pd.read_csv(output / "case_df.tsv", sep="\t")
    case_df = get_subgraph_stats(case_subgraph_dict, "case", bed_files, args["ta_error"])
    case_df.sort_values(["chr", "subgraph", "nodes"]).to_csv(args["output"] / "case_df.tsv", sep="\t", index=False)
    
    control_subgraph_dict = get_subgraphs(args["graph_dir"], "control")
     # control_df = pd.read_csv(output / "control_df.tsv", sep="\t")
    control_df = get_subgraph_stats(control_subgraph_dict, "control", bed_files, args["ta_error"])
    control_df.sort_values(["chr", "subgraph", "nodes"]).to_csv(args["output"] / "control_df.tsv", sep="\t", index=False)
   
    chroms = case_df["chrom"].sort_values().unique()
    
    all_genes_list = []
    # TODO: parallelize this?
    for chrom in chroms:
        # if chrom != "chr6":
        #     continue
        print(chrom)
        
        # get chromosome subsets for case and control
        case_chrom_df = case_df[case_df["chrom"] == chrom]
        control_chrom_df = control_df[control_df["chrom"] == chrom]
        case_chrom_subgraphs = case_subgraph_dict[chrom]
        control_chrom_subgraphs = control_subgraph_dict[chrom]
        annot_chrom_df = annot_df[annot_df["chrom"] == chrom]
        
        # cases as the target
        case_overlaps = pcis_overlaps(case_chrom_df, control_chrom_df)
        case_TA_df, case_overall_df = compare_pcis(case_overlaps, case_chrom_subgraphs, control_chrom_subgraphs)
        case_sig_df = pcis_to_cis(case_overall_df, pval_threshold)
        if len(case_sig_df) != 0:
            case_genes = cis_annotate(case_sig_df, annot_chrom_df)
            case_genes["class"] = "case"
            if verbose > 0:
                # print(len(case_genes))
                print(len(case_genes["marker_symbol"].unique()))
        else:
            case_genes = None
        
        # controls as the target
        control_overlaps = pcis_overlaps(control_chrom_df, case_chrom_df)
        control_TA_df, control_overall_df = compare_pcis(control_overlaps, control_chrom_subgraphs, case_chrom_subgraphs)
        control_sig_df = pcis_to_cis(control_overall_df, pval_threshold)
        if len(control_sig_df) != 0:
            control_genes = cis_annotate(control_sig_df, annot_chrom_df)
            control_genes["class"] = "control"
            if verbose > 0:
                # print(len(control_genes))
                print(len(control_genes["marker_symbol"].unique()))
        else:
            control_genes = None
        
        if case_genes is not None and control_genes is not None:
            both_genes = pd.concat([case_genes, control_genes], ignore_index=True)
        elif case_genes is not None:
            both_genes = case_genes
        elif control_genes is not None:
            both_genes = control_genes
        else:
            pass
        
        if case_genes is not None or control_genes is not None:
            both_genes["chrom"] = chrom
            all_genes_list.append(both_genes)
            if verbose > 0:
                # print(len(both_genes))
                print(len(both_genes["marker_symbol"].unique()))

        # TODO: are there too many repeated genes? Does this make sense that they would be repeated?
        # it appears that sometimes there are multiple CIS in a gene, because the CIS range is quite small
        print(f"""\tsig. genomic features: {both_genes["marker_symbol"].unique().shape[0]}/{annot_chrom_df["Marker Symbol"].unique().shape[0]}""")
        
    # get all genomic features
    all_genes_df = pd.concat(all_genes_list, ignore_index=True)

if __name__ == "__main__": 
    main(load_args())
