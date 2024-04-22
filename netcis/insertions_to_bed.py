import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from docopt import docopt


def load_args():
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


    # # gff3
    # out_df = pd.DataFrame(inserts_df["chr"])
    # out_df.columns = ["seqid"]
    # out_df["source"] = "T2/Onc3"
    # out_df["type"] = "insertion site"
    # out_df["start"] = inserts_df["pos"]
    # out_df["end"] = inserts_df["pos"]
    # out_df["score"] = inserts_df["mapping_quality"]
    # out_df["strand"] = inserts_df["strand"]
    # out_df["phase"] = "."
    # out_df["attributes"] = inserts_df.apply(lambda x: f"ID={x['treatment']}:{x['sampleID']}*{x['sample_sub_id']};color=#{'154360' if x['tpn_promoter_orient'] == '+' else 'D35400'}", axis=1)
    # # strand is still the strand that the insertion was on
    # # tpn_promoter_orient is the orientation w.r.t. the IRL or IRR library
    # # the color in attributes is based on tpn_promoter_orient
    # out_df = pd.concat((out_df, inserts_df.iloc[:, 11:-1]), axis=1)
    # out_df = out_df[out_df["treatment"] == treatment]
    # # out_df.to_csv(insertion_dir.parent / f"{treatment}.gff3", sep="\t", index=False)
    
    
    # # gff2
    # out2_df = pd.DataFrame(inserts_df["chr"])
    # out2_df.columns = ["seqname"]
    # out2_df["source"] = "T2/Onc3"
    # out2_df["feature"] = "insertion site"
    # out2_df["start"] = inserts_df["pos"]
    # out2_df["end"] = inserts_df["pos"]
    # out2_df["score"] = inserts_df["mapping_quality"]
    # out2_df["strand"] = inserts_df["strand"]
    # out2_df["frame"] = "."
    # # out2_df["attributes"] = inserts_df.apply(lambda x: f"ID={x['treatment']}:{x['sampleID']}*{x['sample_sub_id']};color=#{'154360' if x['tpn_promoter_orient'] == '+' else 'D35400'}", axis=1)
    # # out2_df["group "] = inserts_df["sampleID"]
    # out2_df["group "] = inserts_df.apply(lambda x: f"{x['sampleID']}|{x['chr']}|{x['strand']}", axis=1)
    # # strand is still the strand that the insertion was on
    # # tpn_promoter_orient is the orientation w.r.t. the IRL or IRR library
    # # the color in attributes is based on tpn_promoter_orient
    # out2_df["treatment"] = inserts_df["treatment"]
    # out2_df = out2_df[out2_df["treatment"] == treatment]
    # # out2_df.drop("treatment", axis=1).to_csv(insertion_dir.parent / f"{treatment}.gff2", sep="\t", index=False, header=False)
    
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
    out2_df = out2_df.sort_values("chromStart")
    out2_df = out2_df.sort_values("chrom", key=lambda x: x.map(key))
    out2_df.to_csv(insertion_dir.parent / f"{treatment}.bed", sep="\t", index=False, header=False)


if __name__ == "__main__": 
    main(load_args())
    