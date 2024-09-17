import ast, os
from itertools import combinations
from pathlib import Path

from docopt import docopt
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so
from seaborn import axes_style, plotting_context


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
        print(f"number of CIS before: {len(df)}, number of CIS after: {len(data)}")
    return data

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

def remove_Gm(gene_list):
    genes_out = []
    if type(gene_list) is not list:
        gene_list = ast.literal_eval(gene_list)
    for gene in gene_list:
        if not gene.startswith("Gm"):
            genes_out.append(gene)
    return genes_out

def break_cis_rank_ties(df, num_cis):
 
    # values, counts = np.unique(pd.Series([1,2,3,4,4,5,6]).rank(method='min'), return_counts=True)
    # np.unique(counts, return_counts=True)  # if the rankings had no ties, there should be only one value returned
    
    values, counts = np.unique(df["overall_rank"].iloc[:num_cis], return_counts=True)
    ties = values[counts != 1]
    # print(ties)

    df["new_overall_rank"] = df["overall_rank"]
    ranked_cols = ["ranksums_ranked", "fishers_exact_ranked", "total_num_samples_ranked", "total_read_count_ranked"]
    for tie in ties:
        tie_df = df[df["overall_rank"] == tie]
        ind_comb = list(combinations(tie_df.index, 2))
        inc_rank = {}
        for ind1, ind2 in ind_comb:
            ind1_df = tie_df.loc[ind1]
            ind2_df = tie_df.loc[ind2]
            
            a = ind2_df["ranksums_ranked"] - ind1_df["ranksums_ranked"]
            b = ind2_df["fishers_exact_ranked"] - ind1_df["fishers_exact_ranked"]
            c = ind2_df["total_num_samples_ranked"] - ind1_df["total_num_samples_ranked"]
            d = ind2_df["total_read_count_ranked"] - ind1_df["total_read_count_ranked"]
            
            # number of pos or neg numbers. The tiebreak goes to the CIS with more than 2 ranks that are better
            total_rank_diffs = np.sign([a, b, c, d]).sum()
            
            # if the above mthod doesn't work (2 vs. 2) then take the median of each CIS's 4 individual ranks
            ind1_median = tie_df.loc[ind1][ranked_cols].median()
            ind2_median = tie_df.loc[ind2][ranked_cols].median()
            med_res = ind1 if ind1_median > ind2_median else ind2
            
            if total_rank_diffs > 0:  # ind1 is better
                inc_rank[frozenset([ind1, ind2])] = ind2
                if med_res != ind2:
                    # print("interesting 2")  # TODO: could happen that the median does not coincide with whichever CIS is better based on the majority of the 4 rankings
                    pass
            elif total_rank_diffs < 0:  # ind2 is better
                inc_rank[frozenset([ind1, ind2])] = ind1
                if med_res != ind1:
                    # print("interesting 1")  # TODO: could happen that the median does not coincide with whichever CIS is better based on the majority of the 4 rankings
                    pass
            else:  # the difference is 0, use the median
                if ind1_median == ind2_median:
                    # print(f"PANIC, MEDIAN IS THE SAME - {tie}")  # TODO: if median is the same, just take ind1 over ind2 
                    continue
                else:
                    inc_rank[frozenset([ind1, ind2])] = med_res
        
        for key, value in inc_rank.items():
            # for more than 2 samples in a group, each sample is scored against all others
            # therefore, each sample should be incremeted in a way that does not create new ties (not fully tested)
            df.loc[value, 'new_overall_rank'] = df.loc[value, 'new_overall_rank'] + 1
    
    df["new_overall_rank"] = df["new_overall_rank"].astype(int)
    
    values, counts = np.unique(df["new_overall_rank"].iloc[:num_cis], return_counts=True)
    assert len(np.unique(counts)) == 1
    return df

def break_gene_rank_ties(cand_df, annot_df, IS_df):
    # break apart tied genes per CIS by finding which gene start is closest to the most number of insertions
    
    # remove duplicated genes and keep the better scoring CIS gene
    cand_df = cand_df.sort_values(['genes', 'gene_rank'])
    new_cand_genes = cand_df[~cand_df.duplicated('genes')].sort_values('gene_rank').copy()
    
    curr_rank = 1
    value_ties = np.unique(new_cand_genes["overall_rank"])
    for value in value_ties:
        sub_genes_df = new_cand_genes[new_cand_genes["overall_rank"] == value]
        count = len(sub_genes_df)
            
        if count == 1:
            new_cand_genes.loc[sub_genes_df.index, 'gene_rank'] = curr_rank
            curr_rank += 1
            continue
        
        # get all insertion data from case and control 
        case_index_uniq = np.unique(sub_genes_df['case_index'])
        case_index = case_index_uniq[0]
        if type(case_index) is list:
            case_inserts = pd.concat([ IS_df[IS_df['case_index'] == x] for x in case_index ])
        elif pd.isna(case_index):
            case_inserts = pd.DataFrame()
        else:
            case_inserts = IS_df[IS_df['case_index'] == case_index]
        
        control_index_uniq = np.unique(sub_genes_df['control_index'])
        control_index = control_index_uniq[0]
        if type(control_index) is list:
            control_inserts = pd.concat([ IS_df[IS_df['control_index'] == x] for x in control_index ])
        elif pd.isna(control_index):
            control_inserts = pd.DataFrame()
        else:
            control_inserts = IS_df[IS_df['control_index'] == control_index]
        
        sub_genes_inserts = pd.concat([case_inserts, control_inserts])
        
        # determine which gene is closest to the majority of insertions
        res = []
        for i, row in sub_genes_df.iterrows():
            gene_start = annot_df[annot_df['Marker Symbol'] == row['genes']]['genome coordinate start'].to_list()[0]
            # sum of each insertion site's read count (CPM) dividede by the distance from the gene start site
            insert_weight = 1 / np.log((sub_genes_inserts['pos'] - gene_start).abs())
            case_metric = sub_genes_inserts['case_count'] * insert_weight
            control_metric = sub_genes_inserts['control_count'] * insert_weight
            gene_metric = case_metric.sum() + control_metric.sum()
            res.append(pd.DataFrame({'gene': row['genes'], 'metric': gene_metric}, index=[i]))
        gene_rank_df = pd.concat(res).sort_values('metric', ascending=False)
        
        # adjust gene ranks
        for i, row in gene_rank_df.iterrows():
            new_cand_genes.loc[i, 'gene_rank'] = curr_rank
            curr_rank += 1

    # clean up columns to be just the relevant stuff
    cols = [
        'gene_rank', 'genes', "chrom", "CIS_start", "CIS_end", "ranksums", "fishers_exact", 
        "total_num_samples", "case_num_samples", "control_num_samples", "total_read_count",
        "case", "control", "genome_viewer",
        ]
    new_cand_genes = new_cand_genes.sort_values('gene_rank')[cols]
    _, counts = np.unique(new_cand_genes["gene_rank"], return_counts=True)
    assert np.unique(counts)[0] == 1
    return new_cand_genes

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
    
    # TODO: save plot

def edit_tracks_config(track_file):
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
        if header[-4:] == '.gtf':
            
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

def run_gsea():
    # https://gseapy.readthedocs.io/en/latest/introduction.html
    pass


    
def main(args):
    annotation_file: Path = args["annotation"]
    cis_dir: Path = args["CIS_dir"]
    verbose = args["verbose"]
    output: Path = args["output"]
    num_cis = 1000  # TODO: make into another arg

    # load in files
    IS_df = pd.read_csv(cis_dir / "IS.tsv", sep="\t")
    CIS_df = pd.read_csv(cis_dir / "CIS.tsv", sep="\t")
    annot_file_df = pd.read_csv(annotation_file, sep="\t")

    # process results and annotation file
    CIS_df = cast_indexes(CIS_df)
    data_df = filter_low_read_CIS(CIS_df, verbose)
    annot_df = process_annot_file(annot_file_df, args['marker_type'], args['feature_type'], verbose)
    data_df["genes"] = annotate_cis(data_df, annot_df, args['marker_expander'])
    if verbose:
        print(f"number of annotations: {data_df['genes'].apply(lambda x: len(x)).sum()}")
        
    # multi-test correction with Benjaminini-Yekutieli method
    data_df["ranksums_BY"] = false_discovery_control(data_df["ranksums"], method="by")
    data_df["fishers_exact_BY"] = false_discovery_control(data_df["fishers_exact"], method="by")
    data_df["binomial_BY"] = false_discovery_control(data_df["binomial"], method="by")

    # rank CIS by 4 metrics
    data_df["ranksums_ranked"] = data_df["ranksums"].rank()
    data_df["fishers_exact_ranked"] = data_df["fishers_exact"].rank()
    data_df["total_num_samples_ranked"] = data_df["total_num_samples"].rank(ascending=False)
    data_df["total_read_count_ranked"] = data_df["total_read_count"].rank(ascending=False)
    data_df["rank_summed"] = data_df["ranksums_ranked"] + data_df["fishers_exact_ranked"] + data_df["total_num_samples_ranked"] + data_df["total_read_count_ranked"]
    data_df["overall_rank"] = data_df["rank_summed"].rank(method='min')
    data_df = data_df.sort_values("overall_rank")
    data_df = break_cis_rank_ties(data_df, num_cis)

    # add in genome viewer annotation and save CIS data
    data_df["CIS_start"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].min(axis=1).astype(int)
    data_df["CIS_end"] = data_df[["case_pos_min", "case_pos_max", "control_pos_min", "control_pos_max"]].max(axis=1).astype(int)
    data_df["genome_viewer"] = data_df.apply(lambda x: f"""{x["chrom"]}:{x["CIS_start"]}-{x["CIS_end"]}""", axis=1)
    data_df.to_csv(output / "CIS.tsv", sep="\t", index=False)

    # remove G-protein genes cause...there's a lot?
    data_df["genes"] = data_df["genes"].apply(lambda x: remove_Gm(x))

    # generate one gene per row for specified number of CIS
    candidate_df = data_df.iloc[:num_cis].explode("genes").reset_index(drop=True).copy(deep=True)
    candidate_df['gene_rank'] = candidate_df['new_overall_rank'].astype(int)

    # separate results into annotated and unannotated
    not_annotated_df = candidate_df[candidate_df['genes'].isna()]
    annotated_df = candidate_df[~candidate_df['genes'].isna()]

    candidate_genes = break_gene_rank_ties(annotated_df, annot_df, IS_df)
    candidate_genes.to_csv(output / "candidate_genes.tsv", sep="\t", index=False)

    ranked_gene_list = candidate_genes[['genes', 'gene_rank']].copy()
    # GSEA expects genes to be ranked in descending order of best to worst, so this is a simple solution
    ranked_gene_list['gene_rank'] = ranked_gene_list['gene_rank'] * -1
    ranked_gene_list.to_csv(output / "ranked_gene_list.rnk", sep="\t", header=False, index=False)
    ranked_gene_list.drop('gene_rank', axis=1).to_csv(output / "gene_list.tsv", sep="\t", header=False, index=False)
    
    # 4 volcano plots (ranksum, fisher exact, both corrected and uncorrected)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))
    matplot_volcano(ax1, data_df, "ranksums", 0.05, 0, args["case"], args["control"], title="ranksum uncorrected")
    matplot_volcano(ax2, data_df, "ranksums_BY", 0.05, 0, args["case"], args["control"], title="ranksum corrected")
    matplot_volcano(ax3, data_df, "fishers_exact", 0.05, 0, args["case"], args["control"], title="fishers_exact uncorrected")
    matplot_volcano(ax4, data_df, "fishers_exact_BY", 0.05, 0, args["case"], args["control"], title="fishers_exact corrected")
    plt.savefig(output /"volcano_plots.svg")
    
    # pyGenomeViewer genomic tracks
    
    # Gene-set enrichment
    

if __name__ == "__main__":
    main(load_args())
