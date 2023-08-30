import sys, os
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from tqdm import tqdm
from docopt import docopt


def load_args() -> dict:
    doc = """  
    Preprocess .fasta/.fastq files used in Sleeping Beauty transposon mutagenesis screening that used
    both an IRL and IRR library for paired-end reads from Illumina sequencing. This script uses three other
    programs that must be accessible from the command line: cutadapt, bowtie2, and samtools.
    
    Usage: read_depth.py --data DIR --output_prefix STR --bowtie_index DIR --input FILE [options]
    
     -d, --data=DIR                     directory that contains the .fasta or .fastq files of insertions to map. Files can be compressed with gzip
     -o, --output_prefix=DIR            a directory ending with a prefix that will have "-bam" appended to it
     -b, --bowtie_index=DIR             directory that includes the common name of the bowtie2 reference genome index
     -i, --input=FILE                   a file that contains which files to preprocess. See README.md for more information

    Options:
     -h --help                          show this help message and exit
     -v, --verbose=N                    (TODO: NOT USED) print more verbose information using 0, 1 or 2 [default: 0]
     -t, --ntask=INT                    number of threads to use that will speed up cutadapt, bowtie2, and samtools [default: 1]
     -m, --mapq=N                       integer 0-255, higher value is higher quality. See fasta mapQ score for more info [default: 13]
     -n, --npara=INT                    number of parallel processes to preprocess multiple files at the same time [default: 1]
    """

    # remove "--" from args
    new_args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }

    # files and directory args
    new_args["data"] = Path(new_args["data"])
    
    new_args["bam_output"] = Path(new_args["output_prefix"] + "-read_depth_bam")
    new_args["bam_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["bowtie_index"] = Path(new_args["bowtie_index"])
    
    new_args["input"] = Path(new_args["input"])
    
    # int args
    int_opts = ["verbose", "ntask", "npara", "mapq"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])
        
    return new_args

def preprocess_reads(read_f, read_r, mysample_file, ntask, genome_index_dir, mapq_thres) -> None:
    """Process forward and reverse reads: trim transposon and primer, map reads, save to bam file"""
    bam_file = mysample_file.with_suffix(".bam")
    
    mapper = f"bowtie2 -p {ntask} --quiet --local -x {genome_index_dir} -1 {read_f} -2 {read_r}"
    sam_sort = f"samtools sort -@ {ntask} -o {bam_file} 2>/dev/null"
    sam_index = f"samtools index -@ {ntask} {bam_file}"
    os.system(f"{mapper} | {sam_sort}; {sam_index}")
    sys.stdout.flush()
    sys.stderr.flush()

def preprocess_read_helper(iter_args) -> None:
    """helper function for multiprocessing of preprocess_reads()"""
    row, args = iter_args
    
    data_dir = args["data"]
    bam_output_dir = args["bam_output"]
    genome_index_dir = args["bowtie_index"]
    ntask = args["ntask"]
    mapq_thres = args["mapq"]
    
    mysample = row[0]
    irl_F = data_dir / row[1]
    irl_R = data_dir / row[2]
    irl_file = bam_output_dir / (mysample + "_IRL")
    preprocess_reads(irl_F, irl_R, irl_file, ntask, genome_index_dir, mapq_thres)
    
    irr_F = data_dir / row[3]
    irr_R = data_dir / row[4]
    irr_file = bam_output_dir / (mysample + "_IRR")
    preprocess_reads(irr_F, irr_R, irr_file, ntask, genome_index_dir, mapq_thres)

def main() -> None:
    main_args = load_args()
    files_df = pd.read_csv(main_args["input"], sep="\t", header=None)
    iter_args = tqdm([ (row[1], main_args) for row in files_df.iterrows() ])
    with Pool(main_args["npara"]) as p:
        [ x for x in p.imap_unordered(preprocess_read_helper, iter_args) ]
        p.close()
    
if __name__ == "__main__":
    main()
