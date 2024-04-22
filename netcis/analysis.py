import ast
from pathlib import Path

from docopt import docopt
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
import matplotlib.pyplot as plt
import seaborn.objects as so
from seaborn import axes_style


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

def filter_low_read_CIS(df, verbose=0):
    # filter out CIS with low read counts which are highly unlikely to be influential
    # choose threshold that maximizes the removal of 1 nad 2 nubmer of sample CIS, while not removing any with 3 or more
    # in this case that would be 13
    b4_values, b4_counts = np.unique(df["total_num_samples"], return_counts=True)
    cpm_threshold = 0
    old_tmp_str = ""
    total_sample_threshold = 2  # TODO: should this also be a parameter?
    while True:
        af_values, af_counts = np.unique(df[df["total_read_count"] > cpm_threshold]["total_num_samples"], return_counts=True)
        count_differences = 0
        if verbose > 1:
            print(cpm_threshold)
        count_diff_string = ""
        for x in zip(b4_values, b4_counts, af_values, af_counts):
            a, b, c, d = x
            assert a == c
            if b-d != 0:
                count_differences += 1
                tmp_str = f"\tnumber of samples: {a}, CIS filtered out: {b-d}"
                if verbose > 1:
                    print(tmp_str)
                count_diff_string = count_diff_string + tmp_str + '\n'
        if count_differences > total_sample_threshold:
            cpm_threshold -= 1
            break
        cpm_threshold += 1
        old_tmp_str = count_diff_string

    data = df[df["total_read_count"] > cpm_threshold].copy()
    if verbose:
        print(f"at cpm threshold {cpm_threshold}:\n{old_tmp_str}")
        print(f"num CIS before: {len(df)}, num CIS after: {len(data)}")
    return data

def process_annot_file(df, marker_type, feature_type, verbose=0):
    # preprocess annotation file
    # remove genomic features that don't have a genome coordinate start
    df = df[pd.notna(df["genome coordinate start"])]
    # remove unused column
    df = df.drop("Status", axis=1)
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
    
def main(args):
    annotation_file = args["annotation"]
    cis_dir = args["CIS_dir"]
    verbose = args["verbose"]

    # load in files
    # IS_df = pd.read_csv(cis_dir / "IS.tsv", sep="\t")
    CIS_df = pd.read_csv(cis_dir / "CIS.tsv", sep="\t")
    annot_file_df = pd.read_csv(annotation_file, sep="\t")

    # process results and annotation file
    data_df = filter_low_read_CIS(CIS_df, verbose)
    annot_df = process_annot_file(annot_file_df, args['marker_type'], args['feature_type'], verbose)
    data_df["genes"] = annotate_cis(data_df, annot_df, args['marker_expander'])
    if verbose:
        print(f"number of annotations: {data_df['genes'].apply(lambda x: len(x)).sum()}")
        
    # multi-test correction with BY
    data_df["ranksums_BY"] = false_discovery_control(data_df["ranksums"], method="by")
    data_df["fishers_exact_BY"] = false_discovery_control(data_df["fishers_exact"], method="by")
    data_df["binomial_BY"] = false_discovery_control(data_df["binomial"], method="by")

    # rank CIS
    data_df["ranksums_ranked"] = data_df["ranksums"].rank()
    data_df["fishers_exact_ranked"] = data_df["fishers_exact"].rank()
    data_df["total_num_samples_ranked"] = data_df["total_num_samples"].rank(ascending=False)
    data_df["overall_rank"] = (data_df["ranksums_ranked"] + data_df["fishers_exact_ranked"] + data_df["total_num_samples_ranked"]).rank()
    data_df = data_df.sort_values("overall_rank")

    # get candidate annotation list
    data_df["CIS_start"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].min(axis=1)
    data_df["CIS_end"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].max(axis=1)
    data_df["genome_viewer"] = data_df.apply(lambda x: f"""{x["chrom"]}:{x["CIS_start"]}-{x["CIS_end"]}""", axis=1)
    data_df.to_csv(args["output"] / "CIS.tsv", sep="\t", index=False)
        
    # remove G-protein genes cause...there's a lot?
    cols = ["overall_rank", "chrom", "CIS_start", "CIS_end", "genes", "ranksums", "fishers_exact", 
            "total_num_samples", "case_num_samples", "control_num_samples", "case", "control", "genome_viewer"]
    data_df["genes"] = data_df["genes"].apply(lambda x: remove_Gm(x))
    candidate_genes = data_df[cols].explode("genes").reset_index(drop=True)
    candidate_genes.to_csv(args["output"] / "candidate_genes.tsv", sep="\t", index=False)

    # take top genes/features
    top_x_genes = []
    for i, x in enumerate(candidate_genes.itertuples()):
        if not pd.isna(x.genes):
            top_x_genes.append(i)
        if len(top_x_genes) == 200:  # TODO: make another parameter?
            break

    top_genes_df = candidate_genes.iloc[top_x_genes]
    top_genes_df["genes"].to_csv(args["output"] / "top_genes.tsv", sep="\t", index=False, header=False)
    
    # 4 volcano plots (ranksum, fisher exact, both corrected and uncorrected)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))
    matplot_volcano(ax1, data_df, "ranksums", 0.05, 0, args["case"], args["control"], title="ranksum uncorrected")
    matplot_volcano(ax2, data_df, "ranksums_BY", 0.05, 0, args["case"], args["control"], title="ranksum corrected")
    matplot_volcano(ax3, data_df, "fishers_exact", 0.05, 0, args["case"], args["control"], title="fishers_exact uncorrected")
    matplot_volcano(ax4, data_df, "fishers_exact_BY", 0.05, 0, args["case"], args["control"], title="fishers_exact corrected")
    plt.savefig(args["output"] /"volcano_plots.svg")
    
    # genome tracks: convert MGI .rpt file to .bed file
    annot_df["strand"] = annot_df["strand"].fillna(".")
    # BED
    bed_df = pd.DataFrame(annot_df["chrom"])
    bed_df["chromStart"] = annot_df["genome coordinate start"]-1
    bed_df["chromEnd"] = annot_df["genome coordinate end"]
    bed_df["name"] = annot_df["Marker Symbol"]
    bed_df["score"] = 1000
    bed_df["strand"] = annot_df["strand"]
    # bed_df["thickStart"] = annot_df["genome coordinate start"]-1
    # bed_df["thickEnd"] = annot_df["genome coordinate end"]
    # bed_df["itemRGB"] = "255,255,255"
    bed_df.to_csv(args["output"] / "MRK_List2.bed", sep="\t", index=False, header=False)
    

if __name__ == "__main__":
    main(load_args())
