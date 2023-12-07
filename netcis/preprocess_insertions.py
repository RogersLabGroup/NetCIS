from pathlib import Path
import sys
from multiprocessing import Pool

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
     -v, --verbose=N                    print more verbose information using 0, 1 or 2 [default: 0]
     -j, --njobs=N                      an integer for the number of parallel processes to work on multiple files at the same time [default: 1]
    """
    
    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
    
    # files and directory args
    new_args["bam_output"] = Path(new_args["output_prefix"] + "-bam")
    new_args["bam_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["insertions_output"] = Path(new_args["output_prefix"] + "-insertions")
    new_args["insertions_output"].mkdir(parents=True, exist_ok=True)
    new_args["depth_output"] = Path(new_args["output_prefix"] + "-insertions-depth")
    new_args["depth_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["input"] = Path(new_args["input"])
    
    # int args
    new_args["verbose"] = int(new_args["verbose"])
    new_args["njobs"] = int(new_args["njobs"])


    return new_args

def get_insertion_properties(insertion) -> pd.DataFrame:
    """
    record the insertions stats (direction, +/-, and all that)
    NOTE: here is where additional statistics and or properties for each insertion site can be added
    """
    
    if insertion.get_forward_sequence()[:2] == "TA":
        TA_type = "first"
    elif insertion.get_forward_sequence()[-2:] == "TA":
        TA_type = "last"
    else:
        TA_type = "none"
        
    tmp = insertion.get_forward_sequence()
    read_first_last = tmp[:10] + "-" + tmp[-10:]
    tmp = insertion.get_reference_sequence()
    ref_first_last = tmp[:10] + "-" + tmp[-10:]
        
    res = {
        "chr": [insertion.reference_name],
        "pos": [insertion.reference_start + 1],  # 0-based left most coordinate, so we need to add 1 to it to be in line with sam file and gff3 file specs
        "strand": [insertion.is_forward],
        "ref_length": [insertion.reference_length],
        "query_length": [
            insertion.infer_query_length()
        ],  # exludes hard-clipped bases
        "read_length": [
            insertion.infer_read_length()
        ],  # includes hard-clipped bases. should be equal to len(query_sequence)
        "mapping_quality": [insertion.mapping_quality],  # MAPQ: MAPping Quality.
        # MAPQ equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer.
        # A value 255 indicates that the mapping quality is not available.
        # otherwise, the higher the number, the more confident of the quality of the mapping
        # see solution for x in wolfram
        #       254 = -10 * log10(x)
        #       11 = -10 * log10(x)
        "read_name": [insertion.query_name],
        "TA_location": [TA_type],
        "read_first_last": [read_first_last],
        "ref_first_last": [ref_first_last],
        # additional features we could want:
            # query_alignment_sequence - excludes soft clupped bases
            # query_sequence - includes soft clipped bases
    }
    res = pd.DataFrame.from_dict(res)
    return res

def read_is_quality(read) -> bool:
    # that is paired
    if not read.is_paired:
        return False

    # this is mapped
    if not read.is_mapped:
        return False
    
    # # filter reads using quality mapping score
    # if convert_mapq(read.mapping_quality) > mapq_thres:
    #     return False
    
    return True

def process_bam(file, verbose):
    """
    Filter out low quality insertions
    This only can run on paired read sequencing data
    """
    print(file)
    bam = pysam.AlignmentFile(file, "rb")
    # count the coverage for each contig for read normalization later on
    # divide by 2 since these are paired reads and I count two reads that are paired as one read count
    reads_per_chrom_dict = { chrom: int(bam.count(contig=chrom, until_eof=True) / 2) for chrom in bam.references }
    assert bam.unmapped == 0, "there are unmapped reads"
    assert bam.mapped == sum([ bam.count(contig=chrom, until_eof=True) for chrom in bam.references ]), \
        "mapped reads are not equal to total reads, this theoretically shouldn't happen"
    
    insertions = []
    i = 0
    for i, read1 in enumerate(bam.fetch()):  # multiple_iterators=True
        # only look at read 1
        if not read1.is_read1:
            continue
        
        # must have a mate read that is mapped for .mate() to return properly
        if read1.mate_is_unmapped or (not read1.is_paired):
            continue
        
        # read 1 and read 2 must map to the same contig
        read2 = bam.mate(read1)
        if read1.reference_name != read2.reference_name:
            continue
        
        # if the read1 is a quality read, then get the insertions properties
        if read_is_quality(read1):
            insert_properties = get_insertion_properties(read1)
            insertions.append(insert_properties)
                
        # check if read 2 (the mate read) is quality and can be used for insertion properties
        else:  
            if read_is_quality(read2):
                insert_properties = get_insertion_properties(read2)
                insertions.append(insert_properties)
                
    bam.close()
    if verbose > 1:
        print(f"number of reads: {i}")
    
    # check if there were any inseritons at all to avoid errors from pandas.concat()
    if len(insertions) == 0:
        return None, reads_per_chrom_dict
    else:
        df = pd.concat(insertions, axis=0).reset_index(drop=True)
        df["tpn_promoter_orient"] = df["strand"]
        return df, reads_per_chrom_dict

def process_bam_helper(iter_args) -> None:
    mysample, args = iter_args
    bam_dir = args["bam_output"]
    insertions_dir = args["insertions_output"]
    depth_dir = args["depth_output"]
    verbose = args["verbose"]
    
    irl_bam = bam_dir / (mysample + "_IRL.bam")
    irl_sam = bam_dir / (mysample + "_IRL.sam")
    irl_file = pysam.AlignmentFile(irl_sam, "rb")
    irl_total_reads = irl_file.count(until_eof=True) / 2
    irl_file.close()
    # assert irl_total_reads == (irl_file.mapped + irl_file.unmapped)
    
    irr_bam = bam_dir / (mysample + "_IRR.bam")
    irr_sam = bam_dir / (mysample + "_IRR.sam")
    irr_file = pysam.AlignmentFile(irr_sam, "rb")
    irr_total_reads = irr_file.count(until_eof=True) / 2
    irr_file.close()
    # assert irr_total_reads == (irr_file.mapped + irr_file.unmapped)

    # find quality insertion in IRR and IRL libraries and convert them to single insertion site format
    tmp_irl, irl_reads_per_chrom = process_bam(irl_bam, verbose)
    inserts_irl_df = None
    if (tmp_irl is not None):  # if no insertions present, process_bam returns None
        # TODO: normalize read counts (counts per million)
        tmp_irl["count"] = 0  # irrelevant what this holds, it's just a count in the next line
        inserts_irl_df = tmp_irl.groupby(by=["chr", "pos"], sort=False, as_index=False, dropna=False).count()[["chr", "pos", "count"]]
        # divide each insertion site by the total reads on that chromosome
        inserts_irl_df["CPM"] = inserts_irl_df.apply(lambda x: (x["count"] / irl_total_reads) * 1e5, axis=1)
        inserts_irl_df["chrom_norm"] = inserts_irl_df.apply(lambda x: x["count"] / (irl_reads_per_chrom[x["chr"]]), axis=1)
        inserts_irl_df["library"] = "IRL"
        
        # set transposon promoter orientation depending on sequencing library
        # For IRR: + if forward, - if not. For IRL this is reversed. Also, make the orientations easier to read (+/-)
        tmp_irl = tmp_irl.drop("count", axis=1)
        tmp_irl["library"] = "IRL"
        tmp_irl["tpn_promoter_orient"] = ~tmp_irl["tpn_promoter_orient"]
        tmp_irl["strand"] = np.where(tmp_irl["strand"], "+", "-")
        tmp_irl["tpn_promoter_orient"] = np.where(tmp_irl["tpn_promoter_orient"], "+", "-")
        
    tmp_irr, irr_reads_per_chrom = process_bam(irr_bam, verbose)
    inserts_irr_df = None
    if (tmp_irr is not None):
        tmp_irr["count"] = 0
        inserts_irr_df = tmp_irr.groupby(by=["chr", "pos"], sort=False, as_index=False, dropna=False).count()[["chr", "pos", "count"]]
        inserts_irr_df["CPM"] = inserts_irr_df.apply(lambda x: (x["count"] / irr_total_reads) * 1e5, axis=1)
        inserts_irr_df["chrom_norm"] = inserts_irr_df.apply(lambda x: x["count"] / (irr_reads_per_chrom[x["chr"]]), axis=1)
        inserts_irr_df["library"] = "IRR"
        
        tmp_irr = tmp_irr.drop("count", axis=1)
        tmp_irr["library"] = "IRR"
        tmp_irr["strand"] = np.where(tmp_irr["strand"], "+", "-")
        tmp_irr["tpn_promoter_orient"] = np.where(tmp_irr["tpn_promoter_orient"], "+", "-")
    
    # concat of a dataframe and if check if any df is None
    if inserts_irl_df is None and inserts_irr_df is None:
        return
    elif inserts_irl_df is None and inserts_irr_df is not None:
        inserts_df = inserts_irr_df
        individual_inserts = tmp_irr
    elif inserts_irl_df is not None and inserts_irr_df is None:
        individual_inserts = tmp_irl
        inserts_df = inserts_irl_df
    else:
        inserts_df = pd.concat([inserts_irl_df, inserts_irr_df], ignore_index=True)
        inserts_df["TCN"] = inserts_df.apply(lambda x: (x["count"] / (sum(irl_reads_per_chrom.values()) + sum(irr_reads_per_chrom.values()))), axis=1)
        individual_inserts = pd.concat([tmp_irl, tmp_irr], ignore_index=True)


    # # verify that insertions did not count both read1 and read2
    # # do this by checking that the length of 'read names' is the same number as the length of unique read names
    # read_names = inserts_df["read_name"].to_numpy()
    # assert len(np.unique(read_names)) == len(read_names)

    # sort by chr, then pos
    inserts_df = inserts_df.sort_values(["chr", "pos"], ignore_index=True)
    individual_inserts = individual_inserts.sort_values(["chr", "pos"], ignore_index=True)
    
    # TODO: need to update input.tsv with meta info directly below
    # add treatment group and sampleID
    tmp_meta = mysample.split("-")
    if len(tmp_meta) == 3:  # 2020 SB
        inserts_df["treatment"] = tmp_meta[2]
        inserts_df["sampleID"] = tmp_meta[1]
        individual_inserts["treatment"] = tmp_meta[2]
        individual_inserts["sampleID"] = tmp_meta[1]
    elif len(tmp_meta) == 2:  # 2023 SB
        inserts_df["treatment"] = tmp_meta[0]
        inserts_df["sampleID"] = tmp_meta[1]
        individual_inserts["treatment"] = tmp_meta[0]
        individual_inserts["sampleID"] = tmp_meta[1]
    else:  # TODO: gotta change input.tsv to hold extra meta info that I can add
        sys.exit("meta data in mysample is not formmated correctly.")
    
    # save insertions
    inserts_df.to_csv(depth_dir / (mysample + ".tsv"), sep="\t", index=False)
    individual_inserts.to_csv(insertions_dir / (mysample + ".tsv"), sep="\t", index=False)

def main() -> None:
    main_args = load_args()
    files_df = pd.read_csv(main_args["input"], sep="\t", header=None)
    iter_args = tqdm( [(row[1][0], main_args) for row in files_df.iterrows()] )
    with Pool(main_args["njobs"]) as p:
        [ x for x in p.imap_unordered(process_bam_helper, iter_args) ]
        p.close()

if __name__ == "__main__":
    main()
