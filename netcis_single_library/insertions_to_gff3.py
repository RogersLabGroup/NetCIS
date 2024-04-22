import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from docopt import docopt


def load_args():
    """
    Load command-line arguments using docopt.

    Returns:
        dict: Parsed command-line arguments
    """
    
    doc = """  
    Usage: 
        insertions_to_gff3.py --output_prefix DIR --treatment STR [options]
    
     -o, --output_prefix=DIR           a prefix of the output directory that will have "-analysis" appended to it
     -t, --treatment=STR

    Options:
     -h, --help                        show this help message and exit
     -v, --verbose=N                   print more verbose information if available using 0, 1 or 2 [default: 0]

    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
    new_args["insertion_dir"] = Path(new_args["output_prefix"] + "-insertions")
    new_args["verbose"] = int(new_args["verbose"])
    
    if new_args["verbose"] > 1:
        print(new_args)
        
    return new_args
    
def main(args):
    insertion_dir = args["insertion_dir"]
    treatment = args["treatment"]
    
    insert_list = []
    for file in insertion_dir.iterdir():
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

    out_df = pd.DataFrame(inserts_df["chr"])
    out_df.columns = ["seqid"]
    out_df["source"] = "T2/Onc3"
    out_df["type"] = "insertion site"
    out_df["start"] = inserts_df["pos"]
    out_df["end"] = inserts_df["pos"]
    out_df["score"] = inserts_df["mapping_quality"]
    out_df["strand"] = inserts_df["strand"]
    out_df["phase"] = "."
    out_df["attributes"] = inserts_df.apply(lambda x: f"ID={x['treatment']}:{x['sampleID']}*{x['sample_sub_id']};color=#{'154360' if x['tpn_promoter_orient'] == '+' else 'D35400'}", axis=1)
    # strand is still the strand that the insertion was on
    # tpn_promoter_orient is the orientation w.r.t. the IRL or IRR library
    # the color in attributes is based on tpn_promoter_orient
    out_df = pd.concat((out_df, inserts_df.iloc[:, 11:-1]), axis=1)
    out_df = out_df[out_df["treatment"] == treatment]
    out_df.to_csv(insertion_dir.parent / f"{treatment}.gff3", sep="\t", header=False, index=False)

if __name__ == "__main__": 
    main(load_args())
    