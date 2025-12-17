from io import StringIO
from multiprocessing import Pool
from pathlib import Path
import sys, subprocess

import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
from docopt import docopt


def load_args() -> dict:
    doc = """  
    Preprocess the insertions by reading in .bam files of mapped reads. Output files will contain unique
    insertions per row with additional statistics of that insertion. Additionaly, a .tsv file is required
    to idnetify contigs (chromosomes) to keep as well as mapping contigs to a different format if desired.
    
    Usage: preprocess_insertions.py --output_prefix STR --input FILE [options]
    
     -o, --output_prefix=DIR            a directory ending with a prefix that will have "-insertions" appended to it. A directory with "-bam" appended to it should have been created from preprocess_reads.py
     -i, --input=FILE                   a file that contains which files to preprocess, same one used in preprocess_reads.py. See README.md for more information
    
    Options:
     -h, --help                         show this help message and exit
     -v, --verbose=N                    print more verbose information if available using 0, 1 or 2 [default: 0]
     -j, --njobs=N                      an integer for the number of parallel processes to work on multiple files at the same time [default: 1]
    """
    
    # remove "--" from args
    args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
    
    # int args
    args["verbose"] = int(args["verbose"])
    args["njobs"] = int(args["njobs"])
    
    if args["verbose"] > 1:
        print("Arguements given")
        for key, item in args.items():
            print(f"\t{key}: {item} ({type(item)})")
        print("\n")
        
    # files and directory args
    args["bam_output"] = Path(args["output_prefix"] + "-bam")
    args["bam_output"].mkdir(parents=True, exist_ok=True)
    
    args["insertions_dir"] = Path(args["output_prefix"] + "-insertions")
    args["insertions_dir"].mkdir(parents=True, exist_ok=True)
    
    args["strand_dir"] = Path(args["output_prefix"] + "-insertions-strand")
    args["strand_dir"].mkdir(parents=True, exist_ok=True)
    
    args["tpn_orient_dir"] = Path(args["output_prefix"] + "-insertions-strand-tpn_orient")
    args["tpn_orient_dir"].mkdir(parents=True, exist_ok=True)
    
    args["input"] = Path(args["input"])
    
    return args

def sort_chrom_pos(df, chrom, pos):
    """
    Sorts the dataframe by chromosome and position
    """
    key = {
        'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 
        'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 
        'chr11': 11,'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 
        'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 
        'chrX': 20, 'chrY': 21,'chrM': 22,
        }
    custom_sort_key = lambda x: x.map(key)
    df['chrom_custom_sort'] = custom_sort_key(df[chrom])
    df = df.sort_values(by=['chrom_custom_sort', pos], ignore_index=True).drop(columns=['chrom_custom_sort'])
    return df

def run_command(command):
    """
    Run a command in the shell and return the result
    """
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"Error in command: {command}")
        print(result.stderr.decode())
    return result

def process_bam(bam_file):
    """
    Extract insertion information from a BAM file and save it to a TSV file.
    
    *the following code is derived from ChatGPT, verified working - 11/11/2024*
    
    Q:
    how do I extract information from bam file reads? I want the reference name, position, strand, reference length,
    query length, read length, mapping quality, and the read name, the first 10 and last 10 bases of the sequences and
    I only want to use samtools to extract this information.
    
    A: *edited to correct errors*
    samtools view output/2020_SB/results-bam/EL4-18_2-LT_IRL.bam | awk '{
        chr = $3;                                   # Reference name
        pos = $4;                                   # Position
        strand = ($2 == 99) ? "-" : "+";            # Read is paired, proper, mate is revered, and is read1. All reads are read1
        query_length = length($10);                 # Length of the read
        map_quality = $5;                           # Mapping quality
        read_name = $1;                             # Read name
        first_10 = substr($10, 1, 10);              # First 10 bases of read sequence
        last_10 = substr($10, query_length - 9);    # Last 10 bases of read sequence

        print chr, pos, strand, query_length, map_quality, read_name, first_10 "-" last_10 }'

    """
    
    # count the coverage for each contig for read normalization later on
    bam = pysam.AlignmentFile(bam_file, "rb")
    reads_per_chrom_dict = { chrom: int(bam.count(contig=chrom, until_eof=True)) for chrom in bam.references }
    assert bam.unmapped == 0, "there are unmapped reads"
    assert bam.mapped == sum([ bam.count(contig=chrom, until_eof=True) for chrom in bam.references ]), \
        "mapped reads are not equal to total reads, this theoretically shouldn't happen"
    bam.close()
    
    # extract insertion information from bam file
    process_insertion_info = (
        f"""
        samtools view {bam_file} | awk ' OFS="," {{ 
        chr=$3;                                     # Reference name (contig/chromosome)
        pos=$4;                                     # Position
        strand=($2 == 99);                          # Read is paired, proper, mate is reversed, and is read1. All reads are read1 and if read 2 is reversed, read 1 is forward (+)
        query_length=length($10);                   # Length of the read
        map_quality=$5;                             # Mapping quality
        read_name=$1;                               # Read name
        first_10 = substr($10, 1, 10);              # First 10 bases of read sequence
        last_10 = substr($10, query_length - 9);    # Last 10 bases of read sequence
    
        print chr, pos, strand, query_length, map_quality, read_name, first_10 "-" last_10 }}'
        """
    )
    result = run_command(process_insertion_info)
    assert result.returncode == 0, f"Error in process_insertion_info {result.returncode}"
    
    # load data from stdout and add additional features
    cols = ["chr", "pos", "strand", "query_length", "map_quality", "read_name", "read_first_last"]
    df_str = StringIO(result.stdout.decode())
    if len(df_str.getvalue()) == 0:
        return None, reads_per_chrom_dict
    df = pd.read_csv(df_str, sep=",", header=None)
    df.columns = cols
    df['strand'] = df['strand'].astype(bool)
    df["tpn_promoter_orient"] = df["strand"]
    return df, reads_per_chrom_dict
    
def process_individual_insertions(tpn_orient, library, inserts_df):
    if (inserts_df is not None):
        # set transposon promoter orientation depending on sequencing library
        # For IRR: + if forward, - if not. For IRL this is reversed. Also, make the orientations easier to read (+/-)
        if tpn_orient == "-":
            inserts_df["tpn_promoter_orient"] = ~inserts_df["tpn_promoter_orient"]
        inserts_df["tpn_promoter_orient"] = np.where(inserts_df["tpn_promoter_orient"], "+", "-")
        inserts_df["strand"] = np.where(inserts_df["strand"], "+", "-")
        inserts_df["library"] = library
    return inserts_df

def group_insertions(inserts_df, total_reads, reads_per_chrom):
    strand_lib_df = None
    tpn_orient_lib_df = None
    if (inserts_df is not None):
        if len(inserts_df) == 0:
            return strand_lib_df, tpn_orient_lib_df
        # use depth (CPM) for insertion files
        # divide each insertion site by the total reads on that chromosome
        inserts_df["count"] = 0  # irrelevant what this holds, it's just a count in the next line
        strand_lib_df = inserts_df.groupby(
            by=["chr", "pos", "strand", "library"], sort=False, as_index=False, dropna=False
            ).count()[["chr", "pos", "strand", "library", "count"]]
        # TODO: 12/27/24 - should CPM be split by strand?
        strand_lib_df["CPM"] = strand_lib_df.apply(lambda x: (x["count"] / total_reads) * 1e5, axis=1)
        strand_lib_df["chrom_norm"] = strand_lib_df.apply(lambda x: x["count"] / reads_per_chrom[x["chr"]], axis=1)
        
        # use depth (CPM) and strand and transposon orientation for insertion file
        tpn_orient_lib_df = inserts_df.groupby(
            by=["chr", "pos", "strand", "tpn_promoter_orient", "library"], sort=False, as_index=False, dropna=False
            ).count()[["chr", "pos", "strand", "tpn_promoter_orient", "library", "count"]]
        tpn_orient_lib_df["CPM"] = tpn_orient_lib_df.apply(lambda x: (x["count"] / total_reads) * 1e5, axis=1)
        tpn_orient_lib_df["chrom_norm"] = tpn_orient_lib_df.apply(lambda x: x["count"] / reads_per_chrom[x["chr"]], axis=1)
    return strand_lib_df, tpn_orient_lib_df
    
def get_total_reads(prefilter_bam_file):
    # prefiltering bam file has paired reads, while the final bam file has single reads
    file = pysam.AlignmentFile(prefilter_bam_file, "rb")
    total_reads = file.count(until_eof=True) / 2  # when using paired reads, divide by 2
    assert total_reads == (file.mapped + file.unmapped) / 2
    file.close()
    return total_reads

def remove_same_read_name(pos_orient_df, neg_orient_df):
    """remove reads that are in both dataframes"""

    if pos_orient_df is None and neg_orient_df is None:
        return None, None
    elif pos_orient_df is None:
        return None, neg_orient_df
    elif neg_orient_df is None:
        return pos_orient_df, None
    else:
        tmp_pos_orient_df = pos_orient_df[~pos_orient_df['read_name'].isin(neg_orient_df['read_name'])].copy()
        tmp_neg_orient_df = neg_orient_df[~neg_orient_df['read_name'].isin(pos_orient_df['read_name'])].copy()
        
        orient_dup = len(pos_orient_df) - len(tmp_pos_orient_df)
        # if orient_dup:
        #     print(f"\tduplicate reads between pos and neg orientation: {orient_dup}")
            
        return tmp_pos_orient_df, tmp_neg_orient_df


def process_bam_helper(iter_args: dict) -> None:
    """Find quality insertion in IRR and IRL libraries and convert them to single insertion site format.

    Args:
        iter_args (dict): _description_
    """
    
    row, args = iter_args
    sample_id = row.iloc[0]
    bam_dir = args["bam_output"]
    insertions_dir = args["insertions_dir"]
    strand_dir = args["strand_dir"]
    tpn_orient_dir = args["tpn_orient_dir"]
    verbose = args["verbose"]
   
   
    # IRL transposon orientation positive (with strand)
    irl_bam_pos_orient = str(bam_dir / (sample_id + "_IRL.orient_pos.bam"))
    irl_prefilter_pos_orient = str(bam_dir / (sample_id + "_IRL.prefiltering-orient_pos.bam"))
    irl_pos_orient_total_reads = get_total_reads(irl_prefilter_pos_orient)
    tmp_irl_pos_orient, irl_pos_orient_reads_per_chrom = process_bam(irl_bam_pos_orient)
    individual_irl_pos_orient = process_individual_insertions('+', 'IRL', tmp_irl_pos_orient)
    
    # IRL transposon orientation negative (against strand)
    irl_bam_neg_orient = str(bam_dir / (sample_id + "_IRL.orient_neg.bam"))
    irl_prefilter_neg_orient = str(bam_dir / (sample_id + "_IRL.prefiltering-orient_neg.bam"))
    irl_neg_orient_total_reads = get_total_reads(irl_prefilter_neg_orient)
    tmp_irl_neg_orient, irl_neg_orient_reads_per_chrom = process_bam(irl_bam_neg_orient)
    individual_irl_neg_orient = process_individual_insertions('-', 'IRL', tmp_irl_neg_orient)

    # remove same read names and return the new pos and neg orient dataframes for group_insertions
    tmp_ind_irl_pos_orient, tmp_ind_irl_neg_orient = remove_same_read_name(individual_irl_pos_orient, individual_irl_neg_orient)
    strand_irl_df_pos_orient, tpn_orient_irl_df_pos_orient = group_insertions(tmp_ind_irl_pos_orient, irl_pos_orient_total_reads, irl_pos_orient_reads_per_chrom)
    strand_irl_df_neg_orient, tpn_orient_irl_df_neg_orient = group_insertions(tmp_ind_irl_neg_orient, irl_neg_orient_total_reads, irl_neg_orient_reads_per_chrom)


    # IRR transposon orientation positive (with strand)
    irr_bam_pos_orient = str(bam_dir / (sample_id + "_IRR.orient_pos.bam"))
    irr_prefilter_pos_orient = str(bam_dir / (sample_id + "_IRR.prefiltering-orient_pos.bam"))
    irr_pos_orient_total_reads = get_total_reads(irr_prefilter_pos_orient)
    tmp_irr_pos_orient, irr_pos_orient_reads_per_chrom = process_bam(irr_bam_pos_orient)
    individual_irr_pos_orient = process_individual_insertions('+', 'IRR', tmp_irr_pos_orient)
    
    # IRR transposon orientation negative (against strand)
    irr_bam_neg_orient = str(bam_dir / (sample_id + "_IRR.orient_neg.bam"))
    irr_prefilter_neg_orient = str(bam_dir / (sample_id + "_IRR.prefiltering-orient_neg.bam"))
    irr_neg_orient_total_reads = get_total_reads(irr_prefilter_neg_orient)
    tmp_irr_neg_orient, irr_neg_orient_reads_per_chrom = process_bam(irr_bam_neg_orient)
    individual_irr_neg_orient = process_individual_insertions('-', 'IRR', tmp_irr_neg_orient)
    
    # remove same read names and return the new pos and neg orient dataframes for group_insertions
    tmp_ind_irr_pos_orient, tmp_ind_irr_neg_orient = remove_same_read_name(individual_irr_pos_orient, individual_irr_neg_orient)
    strand_irr_df_pos_orient, tpn_orient_irr_df_pos_orient = group_insertions(tmp_ind_irr_pos_orient, irr_pos_orient_total_reads, irr_pos_orient_reads_per_chrom)
    strand_irr_df_neg_orient, tpn_orient_irr_df_neg_orient = group_insertions(tmp_ind_irr_neg_orient, irr_neg_orient_total_reads, irr_neg_orient_reads_per_chrom)
    
    
    # combine library and transposon orientation dataframes
    individual_df = pd.concat([tmp_ind_irl_pos_orient, tmp_ind_irl_neg_orient, tmp_ind_irr_pos_orient, tmp_ind_irr_neg_orient], ignore_index=True)
    strand_df = pd.concat([strand_irl_df_pos_orient, strand_irl_df_neg_orient, strand_irr_df_pos_orient, strand_irr_df_neg_orient], ignore_index=True)
    tpn_orient_df = pd.concat([tpn_orient_irl_df_pos_orient, tpn_orient_irl_df_neg_orient, tpn_orient_irr_df_pos_orient, tpn_orient_irr_df_neg_orient], ignore_index=True)
    
    # sort by chr, then pos
    individual_df = sort_chrom_pos(individual_df, 'chr', 'pos')
    strand_df = sort_chrom_pos(strand_df, 'chr', 'pos')
    tpn_orient_df = sort_chrom_pos(tpn_orient_df, 'chr', 'pos')
    
    # add sample_id, treatment group, and optional metadata
    for i, col in enumerate(row.index):
        if i >= 1 and i <=4:
            continue
        else:
            individual_df[col] = row[col]
            strand_df[col] = row[col]
            tpn_orient_df[col] = row[col]
    
    # save insertions
    individual_df.to_pickle(insertions_dir / (sample_id + ".pkl"))      # for insertions_to_bed.py
    strand_df.to_pickle(strand_dir / (sample_id + ".pkl"))              # for pcis_networks.py
    tpn_orient_df.to_pickle(tpn_orient_dir / (sample_id + ".pkl"))      # for a future version of NetCIS that uses strand and tpn orientation info pCIS


def main() -> None:
    main_args = load_args()
    files_df = pd.read_csv(main_args["input"], sep="\t")
    iter_args = [ (row[1], main_args) for row in files_df.iterrows() ]
    if main_args['verbose']:
        print("preprocess_insertions.py:")
        iter_args = tqdm(iter_args)
        
    with Pool(main_args["njobs"]) as p:
        [ x for x in p.imap_unordered(process_bam_helper, iter_args, chunksize=1) ]
        p.close()
        p.join()

    if main_args['verbose']:
        print()
        
if __name__ == "__main__":
    main()
