import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from docopt import docopt


def load_args():
    doc = """  
    Usage: 
        insertions_to_bed.py --output_prefix DIR --treatment STR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-analysis" appended to it
     -t, --treatment=STR

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
    
    return args
    
def main(args):
    insertion_dir = args["insertion_dir"]
    treatment = args["treatment"]
    # TODO: since we are going to be reading in meta files, this script should work without
    # a given treatment arg and simply run for the treatment given in the meta file
    insert_list = []
    for file in tqdm(insertion_dir.iterdir()):
        tmp_df = pd.read_csv(file, sep="\t")
        tmp_meta = file.name.split(".")[0].split("-")
        # TODO: this needs to be cleaned up and a better standard implemented.
        # maybe a standardized meta-data file?
        if len(tmp_meta) == 3:  # 2020 SB
            tmp_df["treatment"] = tmp_meta[2]
            tmp_df["sampleID"] = tmp_meta[1]
            tmp_df["cell_type"] = tmp_meta[0]
        elif len(tmp_meta) == 2:  # 2023 SB
            tmp_df["treatment"] = tmp_meta[0]
            tmp_df["sampleID"] = tmp_meta[1]
            
        tmp_df["sample_sub_id"] = range(len(tmp_df))
        insert_list.append(tmp_df)
    inserts_df = pd.concat(insert_list, ignore_index=True)

    key = {
        'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 
        'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 
        'chr11': 11,'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 
        'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 
        'chrX': 20, 'chrY': 21,'chrM': 22}
    
    # bed
    out2_df = pd.DataFrame(inserts_df["chr"])
    out2_df.columns = ["chrom"]
    out2_df["chromStart"] = inserts_df["pos"]-1  # 0-index based
    out2_df["chromEnd"] = inserts_df["pos"]  # ending not inclusive
    out2_df["name"] = inserts_df["treatment"]  # inserts_df["sampleID"] or "test"
    out2_df["score"] = 1000
    out2_df["strand"] = inserts_df["strand"]
    out2_df["thickStart"] = out2_df["chromStart"]-1
    out2_df["thickEnd"] = out2_df["chromStart"]
    out2_df["itemRGB"] = "255,0,0"  # TODO: change color for insertion direction
    # out2_df["blockCount"] = ""
    # out2_df["blockSizes"] = ""
    # out2_df["blockStart"] = ""
    out2_df = out2_df[out2_df["name"] == treatment]
    out2_df = out2_df.sort_values("chrom", key=lambda x: x.map(key)).sort_values("chromStart")
    out2_df.to_csv(insertion_dir.parent / f"{treatment}.bed", sep="\t", index=False, header=False)


if __name__ == "__main__": 
    main(load_args())
    