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
    
    Usage: preprocess_reads.py --data DIR --output_prefix STR --bowtie_index DIR --input FILE --irl STR --irr STR --primer STR [options]
    
     -d, --data=DIR                     directory that contains the .fasta or .fastq files of insertions to map. Files can be compressed with gzip
     -o, --output_prefix=DIR            a directory ending with a prefix that will have "-bam" appended to it
     -b, --bowtie_index=DIR             directory that includes the common name of the bowtie2 reference genome index
     -i, --input=FILE                   a file that contains which files to preprocess. See README.md for more information
     -l, --irl=STR                      a string that is the 5'-3' transposon sequence for the IRL library
     -r, --irr=STR                      a string that is the 5'-3' transposon sequence for the IRR library
     -p, --primer=STR                   a string that is the 5'-3' primer sequence

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
    
    new_args["bam_output"] = Path(new_args["output_prefix"] + "-bam")
    new_args["bam_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["report_output"] = new_args["bam_output"].with_name("reports-preprocessing")
    new_args["report_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["bowtie_index"] = Path(new_args["bowtie_index"])
    
    new_args["input"] = Path(new_args["input"])
    
    # sequence args
    new_args["irl"] = Seq(new_args["irl"])
    new_args["irr"] = Seq(new_args["irr"])
    new_args["primer"] = Seq(new_args["primer"])
    
    # int args
    int_opts = ["verbose", "ntask", "npara", "mapq"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])
        
    return new_args

def preprocess_reads(tpn, primer, read_f, read_r, mysample_file, ntask, genome_index_dir, mapq_thres, report_output) -> None:
    """Process forward and reverse reads: trim transposon and primer, map reads, save to bam file"""
    tpn_c = tpn.complement()
    primer_c = primer.complement()
    
    trim1_f = mysample_file.with_name("trim1-" + read_f.name)
    trim1_r = mysample_file.with_name("trim1-" + read_r.name)
    pre_bam_file = mysample_file.with_suffix(".prefiltering.bam")
    bam_file = mysample_file.with_suffix(".bam")
    
    cutadapt_report = mysample_file.with_suffix(".cutadapt.txt").name
    bowtie_report = mysample_file.with_suffix(".bowtie.txt").name
    bam_report = mysample_file.with_suffix(".idxstats.txt").name
    
    # trim reads: for an explanation on trimming Illumina paired-end transposon reads see https://github.com/marcelm/cutadapt/issues/711 
    cutadapt = f"cutadapt -j {ntask} --discard-untrimmed -a {tpn}...{primer} -A {primer_c}...{tpn_c} -o {trim1_f} -p {trim1_r} {read_f} {read_r} > {report_output / cutadapt_report} 2> /dev/null"

    # map reads
    mapper = f"bowtie2 -p {ntask} --very-sensitive-local --local -x {genome_index_dir} -1 {trim1_f} -2 {trim1_r} -S {pre_bam_file} 2> {report_output / bowtie_report}"
    
    # filter reads
    sam_sort = f"samtools sort -@ {ntask} -o {pre_bam_file} {pre_bam_file} > /dev/null 2>&1"
    sam_filter = f"samtools view -h -@ {ntask} -f 2 -q {mapq_thres} -1 -o {bam_file} {pre_bam_file}"
    sam_index = f"samtools index -@ {ntask} {bam_file}"
    sam_stats = f"samtools idxstats {bam_file} > {report_output / bam_report}"
    os.system(f"{cutadapt}; {mapper}; {sam_sort}; {sam_filter}; {sam_index}; {sam_stats}")
    os.system(f"rm {trim1_f} {trim1_r} {pre_bam_file}")
    sys.stdout.flush()
    sys.stderr.flush()

def preprocess_read_helper(iter_args) -> None:
    """helper function for multiprocessing of preprocess_reads()"""
    row, args = iter_args
    
    data_dir = args["data"]
    bam_output_dir = args["bam_output"]
    genome_index_dir = args["bowtie_index"]
    ntask = args["ntask"]
    irl_tpn = args["irl"]
    irr_tpn = args["irr"]
    primer = args["primer"]
    mapq_thres = args["mapq"]
    report_output_dir = args["report_output"]
    
    mysample = row[0]
    irl_F = data_dir / row[1]
    irl_R = data_dir / row[2]
    irl_file = bam_output_dir / (mysample + "_IRL")
    preprocess_reads(irl_tpn, primer, irl_F, irl_R, irl_file, ntask, genome_index_dir, mapq_thres, report_output_dir)
    
    irr_F = data_dir / row[3]
    irr_R = data_dir / row[4]
    irr_file = bam_output_dir / (mysample + "_IRR")
    preprocess_reads(irr_tpn, primer, irr_F, irr_R, irr_file, ntask, genome_index_dir, mapq_thres, report_output_dir)

def main() -> None:
    main_args = load_args()
    files_df = pd.read_csv(main_args["input"], sep="\t", header=None)
    iter_args = tqdm([ (row[1], main_args) for row in files_df.iterrows() ])
    with Pool(main_args["npara"]) as p:
        [ x for x in p.imap_unordered(preprocess_read_helper, iter_args, chunksize=1) ]
        p.close()
    
if __name__ == "__main__":
    main()
