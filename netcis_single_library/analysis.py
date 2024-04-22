from pathlib import Path

from docopt import docopt
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control


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
     -g, --gene_annot=FILE             MGI's mouse menetic markers excluding withdrawn genes
     
    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information if available using 0, 1 or 2 [default: 0]
     -p, --pval_threshold=N            p-value to exclude pCIS for significance [default: 0.05]
     -x, --gene_expander=N             number of base pairs to extend the boundaries of genes when annotation genes to CIS [default: 50000]
     -j, --njobs=N                     number of processes to run [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    int_opts = ["verbose", "njobs", "gene_expander"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])
    
    new_args["pval_threshold"] = float(new_args["pval_threshold"])
    
    new_args["CIS_dir"] = Path(new_args["output_prefix"] + "-CIS")
    new_args["gene_annot"] = Path(new_args["gene_annot"])
    new_args["output"] = Path(new_args["output_prefix"] + "-analysis")
    new_args["output"].mkdir(exist_ok=True)
    
    if new_args["verbose"] > 1:
        print(new_args)
        
    return new_args

def main(args):
    pass



    annot_df = pd.read_csv(args["gene_annot"], sep="\t")
    annot_df = annot_df[pd.notna(annot_df["genome coordinate start"])].drop("Status", axis=1)
    annot_df["chrom"] = annot_df["Chr"].apply(lambda x: f"chr{x}")
    annot_df = annot_df.sort_values("chrom")
    # print(annot_df["Marker Type"].unique())
    # TODO: what about the strand in annot_df?
    


    
    
    # TODO: Annotate pCIS
    # do this per chrom
    # trim down annotation dataframe to just genes
    annot_chrom_genes = annot_chrom_df[annot_chrom_df["Marker Type"] == "Gene"]
    gene_names = annot_chrom_genes["Marker Symbol"].to_numpy()
    # find genes within a pCIS that includes the gene expander range
    pos_min = CIS_df[["case_pos_min", "control_pos_min"]].min(axis=1).to_numpy().reshape(-1, 1)
    pos_max = CIS_df[["case_pos_max", "control_pos_max"]].max(axis=1).to_numpy().reshape(-1, 1)
    gene_start = (annot_chrom_genes["genome coordinate start"] - gene_expander).to_numpy().reshape(1, -1)
    gene_end = (annot_chrom_genes["genome coordinate end"] + gene_expander).to_numpy().reshape(1, -1)
    tmp = (pos_min <= gene_end) & (pos_max >= gene_start)
    # add on genes to pCIS
    CIS_df["genes"] = [ list(gene_names[tmp[i]]) for i in range(tmp.shape[0]) ]


    
    # TODO: multi-test correction with BY
    CIS_df["ranksums_BY"] = false_discovery_control(CIS_df["ranksums"], method="by")
    CIS_df["fishers_exact_BY"] = false_discovery_control(CIS_df["fishers_exact"], method="by")
    CIS_df["binomial_BY"] = false_discovery_control(CIS_df["binomial"], method="by")
    
    # TODO: rank CIS
    

    # TODO: filter down to significant CIS
    # fet = CIS_df["fishers_exact"] <= pval_threshold
    # rst = CIS_df["ranksums"] <= pval_threshold
    # bit = CIS_df["binomial"] <= pval_threshold
    # CIS_df = CIS_df[ fet | rst | bit ]
    
    
    # TODO: get candidate annotation list?
    
    
if __name__ == "__main__": 
    main(load_args())