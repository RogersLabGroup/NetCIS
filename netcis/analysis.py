import ast, os, pickle, sys
from pathlib import Path

from docopt import docopt
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import pairwise_distances
import ranky as rk
import networkx as nx
import gseapy as gp

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import seaborn.objects as so
from seaborn import axes_style, plotting_context


def load_args() -> dict:
    doc = """
    Analyze the CIS results with multiple test corrections, genomic annotations, and plots

    Usage:
        analysis.py --output_prefix DIR --case STR --control STR --annotation FILE --gene_sets FILE [options]

     -o, --output_prefix=DIR            a prefix of the output directory that will have "-analysis" appended to it
     -a, --case=STR                     treatment type value to use as case
     -b, --control=STR                  treatment type value to use as control
     -g, --annotation=FILE              MGI's mouse menetic markers excluding withdrawn genes
     -s, --gene_sets=FILE               Pathway gene sets as a gmt file
     -t, --threshold=N                  maximum distance to connect two insertions together in the CIS network. We suggest not going over the default value [default: 50000]


    Options:
     -h, --help                         show this help message and exit
     -v, --verbose=N                    print more verbose information if available using 0, 1 or 2 [default: 0]
     -p, --pval_threshold=N             a float p-value to exclude pCIS for significance [default: 0.05]
     -x, --marker_expander=N            an integer number of base pairs to extend the boundaries of genomic annotations [default: 5000]
     -m, --marker_type=STR              marker type to annotate based on MRK_List2.rpt file. An empty string "" will not do any filtering [default: Gene]
     -f, --feature_type=STR             marker feature type to annotate based on MRK_List2.rpt file. An empty string "" will not do any filtering [default: protein coding gene]
     -n, --num_cis=N                    and integer number of top ranked CIS to create genome viewer plots for. [default: 20]
     -j, --sim_thresh=N                 theshold for the Jaccard distance used in removing redundant pathways for gene set enrichment [default: 0.5]
    """
    
    # remove "--" from args
    args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
        
    int_opts = ["marker_expander", "verbose", 'num_cis', "threshold"]
    for opts in int_opts:
        args[opts] = int(args[opts])
    
    float_opts = ["pval_threshold", "sim_thresh"]
    for opts in float_opts:
        args[opts] = float(args[opts])
    
    if args["verbose"] > 1:
        print("Arguements given")
        for key, item in args.items():
            print(f"\t{key}: {item}")
        print("\n")
        
    args["CIS_dir"] = Path(args["output_prefix"] + "-CIS") / f"{args['case']}-{args['control']}" / f"{args['threshold']}"
    args["annotation"] = Path(args["annotation"])
    args["output"] = Path(args["output_prefix"] + "-analysis") / f"{args['case']}-{args['control']}" / f"{args['threshold']}"
    args["output"].mkdir(exist_ok=True, parents=True)

    return args

def sort_chrom_pos(df, chrom, pos):
    key = {
        'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 
        'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 
        'chr11': 11,'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 
        'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 
        'chrX': 20, 'chrY': 21,'chrM': 22,
        }
    df = df.copy(deep=True)
    custom_sort_key = lambda x: x.map(key)
    df['chrom_custom_sort'] = custom_sort_key(df[chrom])
    df = df.sort_values(by=['chrom_custom_sort', pos], ignore_index=True).drop(columns=['chrom_custom_sort'])
    return df

def cast_indexes(df):
    # some of the indexes are a string representation of a list and require casting to list type
    for i, row in df.iterrows():
        
        case_ind = row['case_index']
        if type(case_ind) is str:
            df.at[i, 'case_index'] = ast.literal_eval(case_ind)

        control_ind = row['control_index']
        if type(control_ind) is str:
            df.at[i, 'control_index'] = ast.literal_eval(control_ind)
    
    return df

def process_annot_file(df, marker_type, feature_type, verbose=0):
    # preprocess annotation file
    # remove genomic features that don't have a genome coordinate start
    df = df[pd.notna(df["genome coordinate start"])]
    # remove unused column
    df = df.drop("Status", axis=1)
    # remove genes that are heritable phenotypic markers to NaN
    df = df[df['Feature Type'] != 'heritable phenotypic marker']
    # transform Chr column into "chr1" format and sort by Chr
    df["chrom"] = df["Chr"].apply(lambda x: f"chr{x}")
    
    # trim down annotation dataframe to just the chosen marker and feature type
    if verbose:
        print(f"\nmarker types: {df['Marker Type'].unique()}")
    if marker_type != "":
        df = df[(df["Marker Type"] == marker_type)]
    
    if verbose:
        print(f"\nfeature types given marker type is '{marker_type}': {df['Feature Type'].unique()}")
    if feature_type != "":
        df = df[(df["Feature Type"] == feature_type)]

    if verbose:
        print(f"\nnumber of possible annotation: {len(df)}")
    
    df = sort_chrom_pos(df, "chrom", "genome coordinate start")
    df["genome coordinate start"] = df["genome coordinate start"].apply(int)
    df["genome coordinate end"] = df["genome coordinate end"].apply(int)
    # NOTE: for next version
    # what about the strand in annot_df? Does this need to be considered for a CIS by separating by strand?
    return df

def annotate_cis(cis_df, annot_df, marker_expander):
    # Annotate CIS that includes the marker expander range and proper chromosomes
    pos_min = cis_df[["case_pos_min", "control_pos_min"]].min(axis=1).to_numpy().reshape(-1, 1)
    pos_max = cis_df[["case_pos_max", "control_pos_max"]].max(axis=1).to_numpy().reshape(-1, 1)
    cis_chrom = cis_df[["chrom"]].to_numpy().reshape(-1, 1)

    # add marker expander to both ends of the gene irregardless of strand
    marker_start = (annot_df["genome coordinate start"] - marker_expander).to_numpy().reshape(1, -1)
    marker_end = (annot_df["genome coordinate end"] + marker_expander).to_numpy().reshape(1, -1)
    marker_chrom = annot_df["chrom"].to_numpy().reshape(1, -1)

    tmp = (pos_min <= marker_end) & (pos_max >= marker_start) & (cis_chrom == marker_chrom)
    return [ list(annot_df["Marker Symbol"][tmp[i]]) for i in range(tmp.shape[0]) ]

def remove_Gm(gene_list):
    genes_out = []
    if type(gene_list) is not list:
        gene_list = ast.literal_eval(gene_list)
    for gene in gene_list:
        if not gene.startswith("Gm"):
            genes_out.append(gene)
    return genes_out

def volcano_plot(df, pval, pval_thresh, lfc_thresh, title=""):
    data = df.copy(deep=True)
    thres = np.log10(pval_thresh) * -1
    data[pval] = np.log10(data[pval]) * -1
    g = (
        so.Plot(data, x="LFC", y=pval, pointsize=pval)
        # .add(so.Dots(), so.Jitter(1))
        .add(so.Dots(color="grey"), data=data.query(f"{pval} < {thres}"))
        .add(so.Dots(color="blue"), data=data.query(f"{pval} >= {thres}"))
        .scale(y="log")  # type: ignore
        .label(title=title)
    )
    return g

def matplot_volcano(ax, df, pval, pval_thresh, lfc_thresh, case, control, title="", add_text=False):
    # https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html
    ax.scatter(x=df["LFC"], y=df[pval].apply(lambda x:-np.log10(x)), s=3, label="Not significant", color="grey")

    # highlight down- or up- regulated genes
    down = df[ (df['LFC'] <= -lfc_thresh) & (df[pval] <= pval_thresh) ]
    up = df[ (df['LFC'] >= lfc_thresh) & (df[pval] <= pval_thresh) ]

    ax.scatter(x=down['LFC'], y=down[pval].apply(lambda x: -np.log10(x)), s=10, label=f"{control}", color="blue")
    ax.scatter(x=up['LFC'], y=up[pval].apply(lambda x: -np.log10(x)), s=10, label=f"{case}", color="gold")

    if add_text:
        texts_up=[]
        for i,r in up.iterrows():
            texts_up.append(ax.text(x=r['LFC'], y=-np.log10(r[pval]), s=r["genes"]))
        # adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

        texts_down=[]
        for i,r in down.iterrows():
            texts_down.append(ax.text(x=r['LFC'], y=-np.log10(r[pval]), s=r["genes"]))
        # adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    ax.set_xlabel("log2(FC)")
    ax.set_ylabel("-log10(pval)")
    ax.set_title(title)
    if lfc_thresh != 0:
        ax.axvline(-lfc_thresh, color="grey", linestyle="--")
        ax.axvline(lfc_thresh, color="grey", linestyle="--")
    # else:
        # ax.axvline(lfc_thresh, color="grey", linestyle="--")
    ax.axhline(-np.log10(pval_thresh), color="grey", linestyle="--")
    ax.legend(title="Enrichment")

def plot_volcanoes(data_df, pval_threshold, case_group, control_group, output, return_fig=False):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))
    matplot_volcano(ax1, data_df, "ranksums", pval_threshold, 0, case_group, control_group, title="Rank-sum uncorrected")
    matplot_volcano(ax2, data_df, "ranksums-BY", pval_threshold, 0, case_group, control_group, title="Rank-sum corrected")
    matplot_volcano(ax3, data_df, "fishers_exact", pval_threshold, 0, case_group, control_group, title="Fishers exact uncorrected")
    matplot_volcano(ax4, data_df, "fishers_exact-BY", pval_threshold, 0, case_group, control_group, title="Fishers exact corrected")
    ax2.set_ylim(ax1.get_ylim())
    ax4.set_ylim(ax3.get_ylim())
    fig.savefig(output /"volcano_plots.pdf")
    fig.savefig(output /"volcano_plots.svg")
    fig.savefig(output /"volcano_plots.png")
    if return_fig:
        return fig
    else:
        plt.close()
    
def edit_tracks_config(track_file):
    """
        sometimes the variables are commented
        sometimes they need to be uncommented and changed
        or just changed
        this will do all of that
    """
    # why comment out parts if you give it a true or false? just make it false...
    # or just make this into JSON and use that. It's way better for a config file.
    # Looks like the maintainers of pyGenomeTracks will change to something better in version 4.0 (next update)

    vars_to_change = {
        # 'title': None,  # make it shorter
        'labels': 'labels = true',  # change to true
        # 'style': 'UCSC',  # uncomment
        'all_labels_inside': 'all_labels_inside = true',  # uncomment to true
        'labels_in_margin': 'labels_in_margin = true',  # uncomment to true
        # 'merge_transcripts': 'true',  # uncomment to true
        # 'merge_overlapping_exons': 'true',  # uncomment to true
        }

    with open(track_file) as f:
        lines = f.readlines()

    header = None
    out_lines = []
    for line in lines:
        # remove leading and trailing white space 
        line_clean = line.strip()
        new_line = line_clean
        
        # skip empty lines
        if line_clean == "":
            out_lines.append(new_line)
            continue
        
        # check if line is a header - this tells us what config we are editing
        if line_clean[0] == '[':
            header = line_clean.strip('[').strip(']')
            
        # check if line is a variable
        if line_clean.find("=") != -1:
            
            # get just the variable by removing comments and whitespace
            curr_var = line_clean.strip('#').split('=')[0].strip()
            
            # check if the current variable is in the variables to change dictionary
            if curr_var in vars_to_change:
                # print(header, curr_var)
                new_line = vars_to_change[curr_var]
                
                # # for special case of style
                # if curr_var == 'style':
                #     if line_clean.find("UCSC") != -1:
                #         new_line = line_clean.strip('#')
                # else:
                #     new_line = vars_to_change[curr_var]
                    
        # if new_line != line_clean:
        #     print(new_line, line_clean)
            
        out_lines.append(new_line)
        
    with open(track_file, 'w') as f:
        f.write('\n'.join(out_lines))

def edit_pyGV_file_names(pyGT_dir, top_df):
    for row in top_df.itertuples():
        file_name: Path = pyGT_dir / f'test_{row.chrom}-{row.CIS_start}-{row.CIS_end}.png'
        is_annot = 'unannot' if not row.genes else 'annot'
        new_name = pyGT_dir / f'{row.rank:03}-{is_annot}-{row.chrom}:{row.CIS_start}-{row.CIS_end}.png'
        if not file_name.is_file():
            print(f"error: file not found\n\t{file_name}")
        else:
            file_name.rename(new_name)


def get_dist(m, i, j):
    # The metric dist(u=X[i], v=X[j]) is computed and stored in a condensed array whose entry is i < j < m
    # see scipy documentation for pdist for more info
    return m * i + j - ((i + 2) * (i + 1)) // 2

def prepare_gene_set(gene_set_file, gene_set_output, sim_thresh=0.5, verbose=False):
    # it can be time consuming removing redundant pathways so if the gene sets are already made then skip this step
    final_gene_set_name = f'final_gene_sets-sim_{sim_thresh}.pkl'
    if (gene_set_output / final_gene_set_name).is_file():
        if verbose:
            print('gene-set files already exist. Loading...')
        with open(gene_set_output / 'gene_sets.pkl', 'rb') as file1: 
            gene_sets = pickle.load(file1)
        with open(gene_set_output / 'filtered_gene_sets.pkl', 'rb') as file2: 
            filtered_gene_sets = pickle.load(file2)
        with open(gene_set_output / final_gene_set_name, 'rb') as file3: 
            final_gene_sets = pickle.load(file3)
        if verbose:
            print(f"gene-set count: {len(gene_sets)}")
            print(f"filtered gene-set count: {len(filtered_gene_sets)}")
            print(f"final gene-set count: {len(final_gene_sets)}")
        return gene_sets, filtered_gene_sets, final_gene_sets
        
        
    # read gmt file
    gene_set_df = pd.read_csv(gene_set_file, header=None)
    gene_set_df = gene_set_df[0].str.split('\t', expand=True, n=2)
    gene_set_df.columns = ['pathway', "url", "genes"]
    gene_set_df['genes'] = gene_set_df['genes'].str.split('\t')
    gene_set_df['size'] = gene_set_df['genes'].apply(len)
    gene_sets = { x.pathway: x.genes for x in gene_set_df.drop('size', axis=1).itertuples() }
    if verbose:
        print(f"gene-set count: {len(gene_sets)}")
    
    
    # filter pathway list by size of pathway
    min_genes = 10
    max_genes = 300
    filtered_gene_set_df = gene_set_df[(gene_set_df['size'] >= min_genes) & (gene_set_df['size'] <= max_genes)]
    filtered_gene_sets = { x.pathway: x.genes for x in filtered_gene_set_df.drop('size', axis=1).itertuples() }
    if verbose:
        print(f"filtered gene-set count: {len(filtered_gene_sets)}")


    # remove redundant pathways with jaccard distance (dissimilarity) < 0.5
    # the more similar two arrays are the closer to 0 they will be
    # the more distant they are then the closer to 1 they will be


    # make dataframe for gene-set data
    list_genes = sorted(filtered_gene_set_df.explode('genes')['genes'].unique())
    list_pathways = filtered_gene_set_df['pathway'].unique()
    empty_arr = np.zeros((len(list_pathways), len(list_genes)))
    p_g_df = pd.DataFrame(data=empty_arr, index=list_pathways, columns=list_genes)

    # add gene-set data
    for row in filtered_gene_set_df.itertuples():
        p_g_df.loc[row.pathway, row.genes] = 1

    if verbose:
        print('Beginning to remove redundant pathways. This may take a while...', end='')
        sys.stdout.flush()
    # append pathway index if it is the first pathway found to be more similar than the threshold 0.5
    dm2 = pairwise_distances(p_g_df.astype(bool).to_numpy(), metric='jaccard', n_jobs=-1)
    dm = squareform(dm2)
    m = dm2.shape[0]
    removes = []
    for i in range(m):
        for j in range(1, m):
            if i < j:
                ind = get_dist(m, i, j)
                if dm[ind] < sim_thresh:
                    removes.append((j))
    if verbose:
        print('done')
        
    # remove redundant pathways
    remove_rows = p_g_df.index.values[list(set(removes))]
    final_gene_set_df = filtered_gene_set_df[~filtered_gene_set_df['pathway'].isin(remove_rows)]
    final_gene_sets = { x.pathway: x.genes for x in final_gene_set_df.drop('size', axis=1).itertuples() }
    if verbose:
        print(f"final gene-set count: {len(final_gene_sets)}")
    
    with open(gene_set_output / 'gene_sets.pkl', 'wb') as file1: 
        pickle.dump(gene_sets, file1)
    with open(gene_set_output / 'filtered_gene_sets.pkl', 'wb') as file2: 
        pickle.dump(filtered_gene_sets, file2)
    with open(gene_set_output / final_gene_set_name, 'wb') as file3: 
        pickle.dump(final_gene_sets, file3)
        
    return gene_sets, filtered_gene_sets, final_gene_sets

def dot_plot_gse(df, treatment, output, col):
    fig, ax1 = plt.subplots(figsize=(8, 8))
    gp.dotplot(df, column=col, size=6, top_term=10, figsize=(6,8), title = f"Enrichement: {treatment}", cmap="viridis_r", ax=ax1)
    fig.savefig(output / f"enrichr-dotplot-{treatment}.png")
    fig.savefig(output / f"enrichr-dotplot-{treatment}.pdf")
    fig.savefig(output / f"enrichr-dotplot-{treatment}.svg")
    return fig

def enrichment_plot_gse(df, treatment, output, col):
        nodes, edges = gp.enrichment_map(df, column=col, top_term=20)
        G = nx.from_pandas_edgelist(edges, source='src_idx', target='targ_idx', edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

        # Add missing node if there is any
        for node in nodes.index:
            if node not in G.nodes():
                G.add_node(node)
                
        fig, ax = plt.subplots(figsize=(8, 8))

        # init node cooridnates
        # pos=nx.layout.shell_layout(G)
        pos=nx.layout.kamada_kawai_layout(G)

        # draw nodes
        node_size = list(nodes.Hits_ratio*1000)
        node_color = list(nodes['P-value'])
        net_nodes = nx.draw_networkx_nodes(G, pos=pos, cmap='RdYlBu', node_color=node_color, node_size=node_size, ax=ax)
        # make legends
        legend1 = ax.legend(
            *net_nodes.legend_elements("sizes", num=4),
            loc="upper right",
            title="Hits ratio\n* 1000",
            bbox_to_anchor=(1.16, 1),
            borderpad=1,
            labelspacing=2.5,
            )
        ax.add_artist(legend1)
        legend2 = ax.legend(*net_nodes.legend_elements("colors"), loc="lower right", title="P-value", bbox_to_anchor=(1.15, 0))

        # draw node labels
        labels = nodes.Term.to_dict()
        nx.draw_networkx_labels(G, pos=pos, labels=labels, font_size=8, ax=ax)

        # draw edges
        edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
        width = list(map(lambda x: x*10, edge_weight))
        nx.draw_networkx_edges(G, pos=pos, width=width, edge_color='#CDDBD4', ax=ax)
        plt.tight_layout()
        
        fig.savefig(output / f"enrichr-netowrk-{treatment}.png")
        fig.savefig(output / f"enrichr-netowrk-{treatment}.pdf")
        fig.savefig(output / f"enrichr-netowrk-{treatment}.svg")
        
        # save to GraphML format if someone wants to use Cytoscape
        nx.write_graphml(G, output / f'{treatment}.graphml')
        
        return fig

def run_gse(candidate_df, treatment, gene_sets, background, output, return_fig=False):
    
    enriched_df = candidate_df[candidate_df['enriched'] == treatment]
    gene_df1 = enriched_df.explode("genes").reset_index(drop=True).copy(deep=True)
    gene_df2 = gene_df1[~gene_df1['genes'].isna()]
    enriched_genes = set(gene_df2['genes'].to_list())
    gene_list = list(enriched_genes)
    print(f"number of enriched genes in {treatment}: {len(enriched_genes)}")
    
    enrichment = gp.enrich(gene_list=gene_list, gene_sets=gene_sets, background=background, outdir=None)
    
    res_df : pd.DataFrame = enrichment.results.sort_values(['Adjusted P-value', 'P-value'])
    res_df.to_csv(output / f"enrichr-results-{treatment}.tsv", sep='\t')
    
    # dot plot
    try:
        fig1 = dot_plot_gse(res_df, treatment, output, "Adjusted P-value")
        plt.close(fig1)
    except Exception as e:
        print(f"Error in gp.dotplot()\n\t{e}")
        print("\tTrying with un-adjusted p-value")
        try:
            fig1 = dot_plot_gse(res_df, treatment, output, "P-value")
            print("success")
            plt.close(fig1)
        except Exception as e:
            print(f"Error in gp.dotplot()\n\t{e}")
            print("No dot plot can be created")
            fig1 = None
    
    # build graph
    try:
        fig2 = enrichment_plot_gse(res_df, treatment, output, "Adjusted P-value")
        plt.close(fig2)
    except Exception as e:
        print(f"Error in gp.enrichment_map()\n\t{e}")
        print("\tTrying with un-adjusted p-value")
        try:
            fig2 = enrichment_plot_gse(res_df, treatment, output, "P-value")
            print("success")
            plt.close(fig2)
        except Exception as e:
            print(f"Error in gp.enrichment_map()\n\t{e}")
            print("No enrichment plot can be created")
            fig2 = None
    
    if return_fig:
        return (fig1, fig2)
    
def main(args):
    cis_dir: Path = args["CIS_dir"]
    output: Path = args["output"]
    case_group = args['case']
    control_group = args['control']
    annotation_file: Path = args["annotation"]
    gene_set_file = args['gene_sets']
    edge_threshold = args['threshold']

    verbose = args["verbose"]
    pval_threshold = args['pval_threshold']
    marker_expander = args['marker_expander']
    marker_type = args['marker_type']
    feature_type = args['feature_type']
    num_cis = args['num_cis']
    sim_thresh = args['sim_thresh']
    
    if verbose:
        print('analysis.py')
        print(f"\tCase: {case_group}, Control: {control_group}, Edge Threshold: {edge_threshold}")

    # load in files
    # IS_df = pd.read_csv(cis_dir / "IS.tsv", sep="\t")
    CIS_df = pd.read_csv(cis_dir / "CIS.tsv", sep="\t")
    annot_file_df = pd.read_csv(annotation_file, sep="\t")
    annot_df = process_annot_file(annot_file_df, marker_type, feature_type, verbose)

    # process results and annotation fileargs
    data_df = cast_indexes(CIS_df)
    # data_df = filter_low_read_CIS(data_df, verbose)
    data_df["genes"] = annotate_cis(data_df, annot_df, marker_expander)
    # remove G-protein genes cause...there's a lot?
    data_df["genes"] = data_df["genes"].apply(lambda x: remove_Gm(x))
    if verbose:
        print(f"number of annotations: {data_df['genes'].apply(lambda x: len(x)).sum()}")
        
    # multi-test correction with Benjaminini-Yekutieli method and the -log10 transformation
    data_df["ranksums-neglog"] = -np.log10(data_df["ranksums"])
    data_df["fishers_exact-neglog"] = -np.log10(data_df["fishers_exact"])
    with np.errstate(divide = 'ignore'):
        data_df["binomial-neglog"] = -np.log10(data_df["binomial"])
    
    data_df["ranksums-BY"] = false_discovery_control(data_df["ranksums"], method="by")
    data_df["fishers_exact-BY"] = false_discovery_control(data_df["fishers_exact"], method="by")
    data_df["binomial-BY"] = false_discovery_control(data_df["binomial"], method="by")

    data_df["ranksums-BY-neglog"] = -np.log10(data_df["ranksums-BY"])
    data_df["fishers_exact-BY-neglog"] = -np.log10(data_df["fishers_exact-BY"])
    with np.errstate(divide = 'ignore'):
        data_df["binomial-BY-neglog"] = -np.log10(data_df["binomial-BY"])

    # rank CIS by rank aggregation
    rankers = data_df[['ranksums-neglog', 'total_read_count', 'total_num_samples']]
    data_df['rank'] = rk.rank(rk.borda(rankers, method='median'), reverse=True)        
        
    # CIS enrichment
    data_df['enriched'] = ""
    data_df.loc[data_df[(data_df['LFC'] < 0)].index, 'enriched'] = control_group
    data_df.loc[data_df[(data_df['LFC'] > 0)].index, 'enriched'] = case_group

    # add in genome viewer annotation and save CIS data
    data_df["CIS_start"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].min(axis=1).astype(int)
    data_df["CIS_end"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].max(axis=1).astype(int)
    data_df["genome_viewer"] = data_df.apply(lambda x: f"""{x["chrom"]}:{x["CIS_start"] - marker_expander}-{x["CIS_end"] + marker_expander}""", axis=1)
    data_df.to_csv(output / "CIS.tsv", sep="\t", index=False)



    ### make plots for analysis of results
    
    # 4 volcano plots (ranksum, fisher exact, both corrected and uncorrected)
    volcano_output = output / "volcano_plots"
    volcano_output.mkdir(exist_ok=True, parents=True)
    plot_volcanoes(data_df, pval_threshold, case_group, control_group, volcano_output)

    # Gene-set enrichment
    # candidate CIS: select the union of sig. rank sum and fisher exact results for our candidates
    candidate_df = data_df[(data_df['ranksums'] < pval_threshold) | (data_df['fishers_exact'] < pval_threshold)].copy(deep=True)
    if verbose:
        print(f"number of candidate CIS: {len(candidate_df)}")

    # candidate genes
    gene_df1 = candidate_df.explode("genes", ignore_index=True)
    gene_df2 = gene_df1[~gene_df1['genes'].isna()]
    candidate_genes = set(gene_df2['genes'].to_list())
    if verbose:
        print(f"number of candidate genes: {len(candidate_genes)}")

    # background gene set
    background_gene_list = data_df.explode("genes", ignore_index=True)['genes'].dropna().unique().tolist()
    if verbose:
        print(f"number of background genes: {len(background_gene_list)}\n")

    # prepare gene sets
    main_res_dir = output.parent.parent.parent
    gene_set_output = main_res_dir / "processed_gene_sets"
    gene_set_output.mkdir(exist_ok=True, parents=True)
    gene_sets, filtered_gene_sets, final_gene_sets = prepare_gene_set(gene_set_file, gene_set_output, verbose=verbose)
    # NOTE: it is suggested to use the final_gene_sets, as these are pathway that have been filtered for gene size
    # and redundant pathways removed.
    
    # case enrichment
    case_output = output / f"gene_set_enrichment-{case_group}"
    case_output.mkdir(exist_ok=True, parents=True)
    run_gse(candidate_df, case_group, final_gene_sets, background_gene_list, case_output)

    # control enrichment
    control_output = output / f"gene_set_enrichment-{control_group}"
    control_output.mkdir(exist_ok=True, parents=True)
    run_gse(candidate_df, control_group, final_gene_sets, background_gene_list, control_output)


    # pyGenomeViewer genomic tracks
    # genome tracks: convert MGI .rpt file to .bed file
    bed_df = pd.DataFrame(annot_df["chrom"])
    bed_df["chromStart"] = annot_df["genome coordinate start"]-1
    bed_df["chromEnd"] = annot_df["genome coordinate end"]
    bed_df["name"] = annot_df["Marker Symbol"]
    bed_df["score"] = 1000
    bed_df["strand"] = annot_df["strand"].fillna(".")
    bed_df.to_csv(main_res_dir / "MRK_List2.bed", sep="\t", index=False, header=False)


    # NOTE: For next version
    # the number of insertions (CPM) isn't shown, just the insertion site. How should I portray this?
    # maybe with color intensity?

    # make bed file for each region to plot in pyGenomeViewer for candidate genes and for top CIS
    pyGT_helper = output / "pyGT_helper"
    pyGT_helper.mkdir(exist_ok=True, parents=True)
    top_CIS_bed_file = pyGT_helper / "top_CIS.bed"

    # sort bed file
    top_CIS = data_df.sort_values('rank').iloc[:num_cis].copy(deep=True)
    top_CIS = sort_chrom_pos(top_CIS, "chrom", "CIS_start")

    # change CIS start and CIS end to reflect the marker expander
    top_CIS["CIS_start"] = top_CIS["CIS_start"] - marker_expander
    top_CIS["CIS_end"] = top_CIS["CIS_end"] + marker_expander

    # make bed file
    top_CIS_bed = pd.DataFrame(top_CIS["chrom"])
    top_CIS_bed["chromStart"] = top_CIS["CIS_start"]
    top_CIS_bed["chromEnd"] = top_CIS["CIS_end"]
    top_CIS_bed["name"] = top_CIS['genes'].apply(lambda x: " | ".join(x))
    top_CIS_bed["score"] = 1000
    top_CIS_bed["strand"] = "."  # NOTE: add strand specificity in next version
    top_CIS_bed["thickStart"] = top_CIS["CIS_start"]
    top_CIS_bed["thickEnd"] = top_CIS["CIS_end"]
    top_CIS_bed["itemRGB"] = ""  # "255,255,255"
    top_CIS_bed.to_csv(top_CIS_bed_file, sep="\t", index=False, header=False)


    # prepare directories for pyGenomeTracks results
    pyGT_CIS = output / "pyGT_top_CIS"
    pyGT_CIS.mkdir(exist_ok=True, parents=True)

    # get bed case and control bed file and gene annotation file made above for pyGenomeTracks
    bed_files_keep = [f"{case_group}_+", f"{case_group}_-", 
                    f"{control_group}_+", f"{control_group}_-", 
                    annotation_file.stem]
    bed_files_list = bed_files_keep

    for file in main_res_dir.iterdir():
        if file.is_file() and file.stem in bed_files_keep:
            ind = bed_files_keep.index(file.stem)
            bed_files_list[ind] = str(file)
    track_files = " ".join(bed_files_list)

    print()
    print("making genome tracks config file")
    track_out = pyGT_helper / "tracks.ini"
    os.system(f"make_tracks_file --trackFiles {track_files} -o {track_out} > /dev/null 2>&1")
    edit_tracks_config(track_out)

    print('making genome track plots')
    os.system(f"pyGenomeTracks --tracks {track_out} --BED {top_CIS_bed_file} --outFileName {pyGT_CIS / 'test.png'} > /dev/null 2>&1")
    edit_pyGV_file_names(pyGT_CIS, top_CIS)
    print()
    print()

if __name__ == "__main__":
    main(load_args())
