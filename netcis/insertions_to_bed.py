from pathlib import Path

import pandas as pd
from tqdm import tqdm
from docopt import docopt


def load_args():
    doc = """  
    Usage: 
        insertions_to_bed.py --output_prefix DIR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-analysis" appended to it

    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information if available using 0, 1 or 2 [default: 0]
    """

    # remove "--" from args
    args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
    
    args["verbose"] = int(args["verbose"])
    
    if args["verbose"] > 1:
        print("Arguements given")
        for key, item in args.items():
            print(f"\t{key}: {item}")
        print("\n")
    
    args["insertion_dir"] = Path(args["output_prefix"] + "-insertions")
    args["strand_dir"] = Path(args["output_prefix"] + "-insertions-strand")
    args["tpn_orient_dir"] = Path(args["output_prefix"] + "-insertions-strand-tpn_orient")
    
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
    df = df.sort_values(by=['chrom_custom_sort', pos], kind="timsort", ignore_index=True).drop(columns=['chrom_custom_sort'])
    return df


def main(args):
    """Generate individual insertion bed files from insertion directory.

    Args:
        args (dict): input arguments from command line. See load_args() for details.
    """
    insertion_dir = args["insertion_dir"]   # for insertions_to_bed.py
    strand_dir = args["strand_dir"]           # for pcis_networks.py
    tpn_orient_dir = args["tpn_orient_dir"]         # for a future version of NetCIS that uses strand and tpn orientation info pCIS
    verbose = args["verbose"]

    if verbose:
        print("insertions_to_bed.py")
        print("Loading insertions from files...", end="")
    insert_list = [ pd.read_pickle(file) for file in tpn_orient_dir.iterdir()  ]
    if verbose:
        print("done")
    inserts_df = pd.concat(insert_list, ignore_index=True)
    # inserts_df.to_csv(insertion_dir.parent / f"all_insertions.tsv", sep="\t", index=False)
    
    # add strand color to be used in bed file based on if the strand and transposon promoter orientation match
    inserts_df["strand_color"] = ''
    # color blind safe blue if they match
    inserts_df.loc[inserts_df[inserts_df['strand'] == inserts_df['tpn_promoter_orient']].index, 'strand_color'] = "44,123,182"
    # color blind safe red if they don't matched
    inserts_df.loc[inserts_df[inserts_df['strand'] != inserts_df['tpn_promoter_orient']].index, 'strand_color'] = "215,25,28"

    # save bed files for each treatment and strand
    for treatment in inserts_df['treatment'].unique():
        for strand in inserts_df['strand'].unique():
            if verbose:
                print(f"Creating bed file for treatment: {treatment}, and strand: {strand}")
            treatment_df = inserts_df[(inserts_df["treatment"] == treatment) & (inserts_df["strand"] == strand)].copy()
            # treatment_strand_grouped_df = treatment_strand_df.groupby(["chr", "pos"])['CPM'].sum()
            treatment_df['CPM_score'] = (treatment_df['CPM'] / treatment_df['CPM'].max()) * 1000
            treatment_df = sort_chrom_pos(treatment_df, 'chr', 'pos')

            out_df = pd.DataFrame(treatment_df["chr"])
            out_df.columns = ["chrom"]
            out_df["chromStart"] = treatment_df["pos"]-1  # 0-index based
            out_df["chromEnd"] = treatment_df["pos"]  # ending not inclusive
            out_df["name"] = ''
            out_df["score"] = treatment_df['CPM_score']
            out_df["strand"] = treatment_df["strand"]
            out_df["thickStart"] = 0  # out_df["chromStart"]-1
            out_df["thickEnd"] = 0  # out_df["chromStart"]
            out_df["itemRGB"] = treatment_df["strand_color"]
            out_df.to_csv(insertion_dir.parent / f"{treatment}_{strand}.bed", sep="\t", index=False, header=False)
    if verbose:
        print()
        
if __name__ == "__main__": 
    main(load_args())
    