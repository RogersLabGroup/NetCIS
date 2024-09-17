import ast, os
from pathlib import Path

from docopt import docopt
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import pairwise_distances
import ranky as rk

import matplotlib.pyplot as plt
import seaborn.objects as so
from seaborn import axes_style, plotting_context

# pyright: reportArgumentType=false

def load_args() -> dict:
    doc = """
    Analyze the CIS results with multiple test corrections, genomic annotations, and plots

    Usage:
        analysis.py --output_prefix DIR --case STR --control STR --annotation FILE [options]

     -o, --output_prefix=DIR            a prefix of the output directory that will have "-analysis" appended to it
     -a, --case=STR                     treatment type value to use as case
     -b, --control=STR                  treatment type value to use as control
     -g, --annotation=FILE              MGI's mouse menetic markers excluding withdrawn genes

    Options:
     -h, --help                         show this help message and exit
     -v, --verbose=N                    print more verbose information if available using 0, 1 or 2 [default: 0]
     -p, --pval_threshold=N             p-value to exclude pCIS for significance [default: 0.05]
     -x, --marker_expander=N            number of base pairs to extend the boundaries of genomic annotations [default: 5000]
     -m, --marker_type=STR              marker type to annotate based on MRK_List2.rpt file. An empty string "" will not do any filtering [default: Gene]
     -f, --feature_type=STR             marker feature type to annotate based on MRK_List2.rpt file. An empty string "" will not do any filtering [default: protein coding gene]
    """

    # remove "--" from args
    args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
        
    int_opts = ["marker_expander", "verbose"]
    for opts in int_opts:
        args[opts] = int(args[opts])
        
    if args["verbose"] > 1:
        print("Arguements given")
        for key, item in args.items():
            print(f"\t{key}: {item}")
        print("\n")
    
    args["pval_threshold"] = float(args["pval_threshold"])
    
    args["CIS_dir"] = Path(args["output_prefix"] + "-CIS") / f"{args['case']}-{args['control']}"
    args["annotation"] = Path(args["annotation"])
    args["output"] = Path(args["output_prefix"] + "-analysis") / f"{args['case']}-{args['control']}"
    args["output"].mkdir(exist_ok=True, parents=True)

    return args

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
    # TODO: need a more flexible way to choose marker and feature types with exclusion and inclusion params
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
        
    df = df.sort_values(["Chr", "genome coordinate start"])
    df["genome coordinate start"] = df["genome coordinate start"].apply(int)
    df["genome coordinate end"] = df["genome coordinate end"].apply(int)
    # TODO: what about the strand in annot_df? Does this need to be considered for a CIS by separating by strand?
    return df

def annotate_cis(cis_df, annot_df, marker_expander):
    # Annotate CIS that includes the marker expander range and proper chromosomes
    pos_min = cis_df[["case_pos_min", "control_pos_min"]].min(axis=1).to_numpy().reshape(-1, 1)
    pos_max = cis_df[["case_pos_max", "control_pos_max"]].max(axis=1).to_numpy().reshape(-1, 1)
    cis_chrom = cis_df[["chrom"]].to_numpy().reshape(-1, 1)

    marker_start = (annot_df["genome coordinate start"] - marker_expander).to_numpy().reshape(1, -1)
    marker_end = (annot_df["genome coordinate end"] + marker_expander).to_numpy().reshape(1, -1)
    marker_chrom = annot_df["chrom"].to_numpy().reshape(1, -1)

    tmp = (pos_min <= marker_end) & (pos_max >= marker_start) & (cis_chrom == marker_chrom)
    return [ list(annot_df["Marker Symbol"][tmp[i]]) for i in range(tmp.shape[0]) ]

def gene_search(df, gene_search):
    return df[df["genes"] == gene_search]

def remove_Gm(gene_list):
    genes_out = []
    if type(gene_list) is not list:
        gene_list = ast.literal_eval(gene_list)
    for gene in gene_list:
        if not gene.startswith("Gm"):
            genes_out.append(gene)
    return genes_out

def quantile_selection(df, column, quant):
    quant_rank = df[column].quantile(quant)
    return df[df[column] < quant_rank].sort_values(column)

def volcano_plot(df, pval, pval_thresh, lfc_thresh, title=""):
    data = df.copy(deep=True)
    thres = np.log10(pval_thresh) * -1
    data[pval] = np.log10(data[pval]) * -1
    g = (
        so.Plot(data, x="LFC", y=pval, pointsize=pval)
        # .add(so.Dots(), so.Jitter(1))
        .add(so.Dots(color="grey"), data=data.query(f"{pval} < {thres}"))
        .add(so.Dots(color="blue"), data=data.query(f"{pval} >= {thres}"))
        .scale(y="log")
        .label(title=title)
    )
    return g

def matplot_volcano(ax, df, pval, pval_thresh, lfc_thresh, case, control, title="", add_text=False):
    # https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html
    ax.scatter(x=df["LFC"], y=df[pval].apply(lambda x:-np.log10(x)), s=3, label="Not significant", color="grey")

    # highlight down- or up- regulated genes
    down = df[ (df['LFC'] <= -lfc_thresh) & (df[pval] <= pval_thresh) ]
    up = df[ (df['LFC'] >= lfc_thresh) & (df[pval] <= pval_thresh) ]

    ax.scatter(x=down['LFC'], y=down[pval].apply(lambda x: -np.log10(x)), s=10, label=f"{case} < {control}", color="blue")
    ax.scatter(x=up['LFC'], y=up[pval].apply(lambda x: -np.log10(x)), s=10, label=f"{case} > {control}", color="gold")

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
    ax.legend()

def edit_tracks_config(track_file):
    # TODO: figure out the final version of the tracks.ini file that I need to make it look good
    """
        sometimes the variables are commented
        sometimes they need to be uncommented and changed
        or just changed
        this will do all of that
    """
    # why comment out parts if you give it a true or false? just make it false...
    # or just make this into JSON and use that. It's way better for a config file.
    # might want to make an issue on the github repo to change this?

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
            
        # if we are in the annotation configs
        if (header is not None) and (header[-4:] == '.gtf'):
            
            # check if line is a variable
            if line_clean.find("=") != -1:
                
                # get just the variable by removing comments and whitespace
                curr_var = line_clean.strip('#').split('=')[0].strip()
                
                # check if the current variable is in the variables to change dictionary
                if curr_var in vars_to_change:
                    
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

def edit_pyGV_file_names(pyGV_dir, top_df):
    for row in top_df.itertuples():
        file_name = pyGV_dir / f'test_{row.chrom}-{row.CIS_start}-{row.CIS_end}.png'
        is_annot = 'unannot' if not row.genes else 'annot'
        new_name = pyGV_dir / f'{row.new_overall_rank:04}-{is_annot}-{row.chrom}-{row.CIS_start}-{row.CIS_end}.png'
        if not file_name.is_file():
            print(f"error: file not found {pyGV_dir / file_name}")
        else:
            file_name.rename(new_name)
            
def main(args):
    cis_dir: Path = args["CIS_dir"]
    output: Path = args["output"]
    case_group = args['case']
    control_group = args['control']
    annotation_file: Path = args["annotation"]

    verbose = args["verbose"]
    pval_threshold = args['pval_threshold']
    marker_expander = args['marker_expander']
    marker_type = args['marker_type']
    feature_type = args['feature_type']
    num_cis = 100  # TODO: add to args
    

    # load in files
    IS_df = pd.read_csv(cis_dir / "IS.tsv", sep="\t")
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
    data_df["binomial-neglog"] = -np.log10(data_df["binomial"])

    data_df["ranksums-BY"] = false_discovery_control(data_df["ranksums"], method="by")
    data_df["fishers_exact-BY"] = false_discovery_control(data_df["fishers_exact"], method="by")
    data_df["binomial-BY"] = false_discovery_control(data_df["binomial"], method="by")

    data_df["ranksums-BY-neglog"] = -np.log10(data_df["ranksums-BY"])
    data_df["fishers_exact-BY-neglog"] = -np.log10(data_df["fishers_exact-BY"])
    data_df["binomial-BY-neglog"] = -np.log10(data_df["binomial-BY"])

    # rank CIS by rank aggregation
    rankers = data_df[['ranksums-neglog', 'fishers_exact-neglog', 'total_read_count', 'total_num_samples']]
    data_df['rank'] = rk.rank(rk.borda(rankers, method='median'), reverse=True)        
        
    # CIS enrichment
    data_df['enriched'] = ""
    data_df.loc[data_df[(data_df['LFC'] < 0)].index, 'enriched'] = case_group
    data_df.loc[data_df[(data_df['LFC'] > 0)].index, 'enriched'] = control_group

    # add in genome viewer annotation and save CIS data
    data_df["CIS_start"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].min(axis=1).astype(int)
    data_df["CIS_end"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].max(axis=1).astype(int)
    data_df["genome_viewer"] = data_df.apply(lambda x: f"""{x["chrom"]}:{x["CIS_start"] - marker_expander}-{x["CIS_end"] + marker_expander}""", axis=1)
    data_df.to_csv(output / "CIS.tsv", sep="\t", index=False)


    # 4 volcano plots (ranksum, fisher exact, both corrected and uncorrected)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))
    matplot_volcano(ax1, data_df, "ranksums", pval_threshold, 0, args["case"], args["control"], title="Rank-sum uncorrected")
    matplot_volcano(ax2, data_df, "ranksums-BY", pval_threshold, 0, args["case"], args["control"], title="Rank-sum corrected")
    matplot_volcano(ax3, data_df, "fishers_exact", pval_threshold, 0, args["case"], args["control"], title="Fishers exact uncorrected")
    matplot_volcano(ax4, data_df, "fishers_exact-BY", pval_threshold, 0, args["case"], args["control"], title="Fishers exact corrected")
    ax2.set_ylim(ax1.get_ylim())
    ax4.set_ylim(ax3.get_ylim())
    fig.savefig(output /"volcano_plots.pdf")
    fig.savefig(output /"volcano_plots.svg")
    fig.savefig(output /"volcano_plots.png")




    ######## pyGenomeViewer genomic tracks
    
    # TODO: use this for pyGenomeViewer?
    # genome tracks: convert MGI .rpt file to .bed file
    # BED
    bed_df = pd.DataFrame(annot_df["chrom"])
    bed_df["chromStart"] = annot_df["genome coordinate start"]-1
    bed_df["chromEnd"] = annot_df["genome coordinate end"]
    bed_df["name"] = annot_df["Marker Symbol"]
    bed_df["score"] = 1000
    bed_df["strand"] = annot_df["strand"].fillna(".")
    # bed_df["thickStart"] = annot_df["genome coordinate start"]-1
    # bed_df["thickEnd"] = annot_df["genome coordinate end"]
    # bed_df["itemRGB"] = "255,255,255"
    bed_df.to_csv(output/ "MRK_List2.bed", sep="\t", index=False, header=False)
    # TODO: move back into 2020_SB output directory?



    # TODO: the number of insertions isn't shown, just the insertion site. How should I portray this?
    # make bed file for each region to plot in pyGenomeViewer for candidate genes and for top CIS
    top_CIS_bed_file = output / "top_CIS-pyGV.bed"
    top_CIS = data_df.sort_values('rank').iloc[:num_cis].copy()
    top_CIS["CIS_start"] = top_CIS["CIS_start"] - marker_expander
    top_CIS["CIS_end"] = top_CIS["CIS_end"] + marker_expander

    top_CIS_bed = pd.DataFrame(top_CIS["chrom"])
    top_CIS_bed["chromStart"] = top_CIS["CIS_start"]
    top_CIS_bed["chromEnd"] = top_CIS["CIS_end"]
    top_CIS_bed["name"] = top_CIS['genes'].apply(lambda x: "_".join(x))  # TODO: change this to CPM?
    top_CIS_bed["score"] = 1000
    top_CIS_bed["strand"] = "."  # TODO: add strand specificity?
    top_CIS_bed.to_csv(top_CIS_bed_file, sep="\t", index=False, header=False)

    # prepare directories for pyGenomeViewer results
    pyGV_CIS = output / "pyGV_top_CIS"
    pyGV_CIS.mkdir(exist_ok=True, parents=True)

    # get bed files programatically - file annotation.gtf.gz should be in the same high level result directory
    # TODO: maybe make a bed file from the rpt file and use that? Or do they supply one already?
    bed_files_list = []
    for file in output.parent.parent.iterdir():
        if file.is_file():
            bed_files_list.append(str(file))
            
    bed_files_list = [
        'output/2020_SB/LT.bed',
        'output/2020_SB/RT.bed',
        'output/2020_SB/S.bed',
        'output/2020_SB/MRK_List2.bed',
        # 'output/2020_SB/gencode.vM35.annotation.gtf.gz',
        ]
    track_files = " ".join(bed_files_list)

    track_out = output / "tracks.ini"

    # TODO: need a better way to make the tracks config file
    # print("make genome track config file")
    # os.system(f"make_tracks_file --trackFiles {track_files} -o {track_out} > /dev/null")

    # print("edit config file")
    # edit_tracks_config(track_out)


    # make genomic track images
    print('make genome track plots')
    os.system(f"pyGenomeTracks --tracks {track_out} --BED {top_CIS_bed_file} --outFileName {pyGV_CIS / 'test.png'} > /dev/null 2>&1")
    # run2 = f"pyGenomeTracks --tracks {track_out} --BED {top_genes_bed_file} --outFileName {pyGV_genes / 'test.png'} > /dev/null 2>&1"
    # run3 = f"pyGenomeTracks --tracks {track_out} --BED {top_unannot_bed_file} --outFileName {pyGV_unannot / 'test.png'} > /dev/null 2>&1"
    # processes = [ subprocess.Popen(program, shell=True) for program in [run1, run2, run3] ]
    # for process in processes:
    #     process.wait()


    # go through png files and rename to relative ranks or genes
    print('editing genome track file names')
    edit_pyGV_file_names(pyGV_CIS, top_CIS)






    ######## Gene-set enrichment
    # select the union of sig. rank sum and fisher exact results for our candidates
    # candidate_df = data_df[(data_df['ranksums'] < 0.05)]
    # candidate_df = data_df[(data_df['fishers_exact'] < 0.05)]
    candidate_df = data_df[(data_df['ranksums'] < pval_threshold) | (data_df['fishers_exact'] < pval_threshold)]
    gene_df1 = candidate_df.explode("genes").reset_index(drop=True).copy(deep=True)
    gene_df2 = gene_df1[~gene_df1['genes'].isna()]
    candidate_genes = set(gene_df2['genes'].to_list())
    print(len(candidate_df))
    print(len(candidate_genes))
    # adding in the fisher's exact adds 4 new genes (1008 to 1012), but one of them is the experimentally validated Sprr1b


    ##### enrichr #####
    import gseapy as gp
    # over-representation analysis using hypergeometric test
    # gp.enrich is local while gp.enrichr is using Enrichr web services


    background_gene_list = data_df['genes'].explode().dropna().unique().tolist()
    print(f"background gene count: {len(background_gene_list)}")

    # TODO: download MSigDB gmt files for analysis
    # read gmt file
    # gene_set = str(output.parent / 'msigdb.v2023.2.Mm.symbols.gmt')
    # gene_set = str(output.parent / 'm2.all.v2023.2.Mm.symbols.gmt')
    # gene_set = str(output.parent / 'm2.cp.v2023.2.Mm.symbols.gmt')
    gene_set = str(output.parent / 'm5.all.v2023.2.Mm.symbols.gmt')
    # gene_set = str(output.parent / 'm5.go.v2023.2.Mm.symbols.gmt')

    gene_set_df = pd.read_csv(gene_set, header=None)
    gene_set_df = gene_set_df[0].str.split('\t', expand=True, n=2)
    gene_set_df.columns = ['pathway', "url", "genes"]
    gene_set_df['genes'] = gene_set_df['genes'].str.split('\t')
    gene_set_df['size'] = gene_set_df['genes'].apply(len)

    gene_set_list = { x.pathway: x.genes for x in gene_set_df.drop('size', axis=1).itertuples() }
    print(f"gene-set count: {len(gene_set_list)}")

    min_genes = 10
    max_genes = 300
    filtered_gene_set_df = gene_set_df[(gene_set_df['size'] >= min_genes) & (gene_set_df['size'] <= max_genes)]
    filtered_gene_set_list = { x.pathway: x.genes for x in filtered_gene_set_df.drop('size', axis=1).itertuples() }
    print(f"filtered gene-set count: {len(filtered_gene_set_list)}")




    # remove redundant pathways with jaccard distance (dissimilarity) < 0.5
    # the more similar two arrays are the closer to 0 they will be
    # the more distant they are then the closer to 1 they will be

    def get_dist(m, i, j):
        # The metric dist(u=X[i], v=X[j]) is computed and stored in a condensed array whose entry is i < j < m
        # see scipy documentation for pdist for more info
        return m * i + j - ((i + 2) * (i + 1)) // 2

    sim_thres = 0.5

    # make dataframe for gene-set data
    list_genes = sorted(filtered_gene_set_df.explode('genes')['genes'].unique())
    list_pathways = filtered_gene_set_df['pathway'].unique()
    empty_arr = np.zeros((len(list_pathways), len(list_genes)))
    p_g_df = pd.DataFrame(data=empty_arr, index=list_pathways, columns=list_genes)

    # add gene-set data
    for row in filtered_gene_set_df.itertuples():
        p_g_df.loc[row.pathway, row.genes] = 1

    # append pathway index if it is the first pathway found to be more similar than the threshold 0.5
    dm2 = pairwise_distances(p_g_df.astype(bool).to_numpy(), metric='jaccard', n_jobs=-1)
    dm = squareform(dm2)
    m = dm2.shape[0]
    removes = []
    for i in range(m):
        for j in range(1, m):
            if i < j:
                ind = get_dist(m, i, j)
                if dm[ind] < sim_thres:
                    removes.append((j))
                    
    # remove redundant pathways
    remove_rows = p_g_df.index.values[list(set(removes))]
    final_gene_set_df = filtered_gene_set_df[~filtered_gene_set_df['pathway'].isin(remove_rows)]
    final_gene_set_list = { x.pathway: x.genes for x in final_gene_set_df.drop('size', axis=1).itertuples() }
    print(f"final gene-set count: {len(final_gene_set_list)}")

    lt_candidate_df = candidate_df[candidate_df['enriched'] == 'LT']
    gene_df1 = lt_candidate_df.explode("genes").reset_index(drop=True).copy(deep=True)
    gene_df2 = gene_df1[~gene_df1['genes'].isna()]
    lt_candidate_genes = set(gene_df2['genes'].to_list())
    print(len(lt_candidate_genes))

    enr_lt = gp.enrich(
        gene_list=list(lt_candidate_genes),
        gene_sets=final_gene_set_list,
        background=background_gene_list,
        outdir=None,
        verbose=1,
        )
    enr_lt.results.sort_values(['Adjusted P-value', 'P-value']).head(20)


    # enr.results.Term = enr.results.Term.str.split(" \(GO").str[0]

    ax = gp.dotplot(enr_lt.results,
                    column="P-value",  # TODO: can't use "Adjusted P-value"
                    size=6,
                    top_term=10,
                    figsize=(6,8),
                    title = "Enrichement: LT",
                    cmap="viridis_r"
                    )
    plt.show()

    import networkx as nx

    # build graph
    nodes, edges = gp.enrichment_map(enr_lt.results,
                                    column="P-value",  # TODO: can't use "Adjusted P-value"
                                    top_term=20,
                                    )
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
    node_size = list(nodes.Hits_ratio * 1000)
    node_color = list(nodes['P-value'])
    nx.draw_networkx_nodes(G, pos=pos, cmap='RdYlBu', node_color=node_color, node_size=node_size)

    # draw node labels
    labels = nodes.Term.to_dict()
    nx.draw_networkx_labels(G, pos=pos, labels=labels, font_size=8)

    # draw edges
    edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
    width = list(map(lambda x: x*10, edge_weight))
    nx.draw_networkx_edges(G, pos=pos, width=width, edge_color='#CDDBD4')

    plt.show()


    # save to GraphML format and upload to Cytoscape
    nx.write_graphml(G, 'path')




    s_candidate_df = candidate_df[candidate_df['enriched'] == 'S']

    gene_df1 = s_candidate_df.explode("genes").reset_index(drop=True).copy(deep=True)
    gene_df2 = gene_df1[~gene_df1['genes'].isna()]
    s_candidate_genes = set(gene_df2['genes'].to_list())
    print(len(s_candidate_genes))

    enr_s = gp.enrich(
        gene_list=list(s_candidate_genes),
        gene_sets=final_gene_set_list,  # filtered_gene_set_list gene_set_list
        background=background_gene_list,
        outdir=None,
        verbose=1,
        )
    enr_s.results.sort_values(['Adjusted P-value', 'P-value']).head(20)


    # enr.results.Term = enr.results.Term.str.split(" \(GO").str[0]

    ax = gp.dotplot(enr_s.results,
                    column="P-value",  # TODO: can't use "Adjusted P-value"
                    size=6,
                    top_term=10,
                    figsize=(6,8),
                    title = "Enrichement: S",
                    cmap="viridis_r"
                    )
    plt.show()

    import networkx as nx

    # build graph
    nodes, edges = gp.enrichment_map(enr_s.results,
                                    column="P-value",  # TODO: can't use "Adjusted P-value"
                                    top_term=20,
                                    )
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
    node_size = list(nodes.Hits_ratio * 1000)
    node_color = list(nodes['P-value'])
    nx.draw_networkx_nodes(G, pos=pos, cmap='RdYlBu', node_color=node_color, node_size=node_size)

    # draw node labels
    labels = nodes.Term.to_dict()
    nx.draw_networkx_labels(G, pos=pos, labels=labels, font_size=8)

    # draw edges
    edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
    width = list(map(lambda x: x*10, edge_weight))
    nx.draw_networkx_edges(G, pos=pos, width=width, edge_color='#CDDBD4')

    plt.show()


    # save to GraphML format and upload to Cytoscape
    nx.write_graphml(G, 'path')



if __name__ == "__main__":
    main(load_args())
