from pathlib import Path
import sys
import os
from multiprocessing import Pool

import pandas as pd
from Bio.Seq import Seq
from tqdm import tqdm
from docopt import docopt


def load_args() -> dict:
    doc = """  
    Preprocess .fasta/.fastq files used in Sleeping Beauty transposon mutagenesis screening that used
    both an IRL and IRR library for paired-end reads from Illumina sequencing. This script uses three other
    programs that must be accessible from the command line: cutadapt, bowtie2, and samtools.
    
    Usage: preprocess_reads.py --data DIR --output_prefix STR --bowtie_index DIR --input FILE --irl_tpn STR --irr_tpn STR --adapter STR [options]
    
     -d, --data=DIR                     directory that contains the .fasta or .fastq files of insertions to map. Files can be compressed with gzip
     -o, --output_prefix=DIR            a directory ending with a prefix that will have "-bam" appended to it
     -b, --bowtie_index=DIR             directory that includes the common name of the bowtie2 reference genome index
     -i, --input=FILE                   a file that contains which files to preprocess. See README.md for more information
     -l, --irl=STR                      a string that is the 5'-3' transposon sequence for the IRL library
     -r, --irr=STR                      a string that is the 5'-3' transposon sequence for the IRR library
     -a, --adapter=STR                  a string that is the 5'-3' adapter sequence

    Options:
     -h --help                      show this help message and exit
     -v, --verbose=N                    (NOT USED) print more verbose information using 0, 1 or 2 [default: 0]
     -t, --ntask=INT                    number of threads to use that will speed up cutadapt, bowtie2, and samtools [default: 1]
     -p, --npara=INT                    number of parallel processes to preprocess multiple files at the same time [default: 1]
    """

    # remove "--" from args
    new_args = { new_args[key.split("-")[-1]]: value for key, value in docopt(doc).items() }

    # files and directory args
    new_args["data"] = Path(new_args["data"])
    
    new_args["bam_output"] = Path(new_args["output_prefix"] + "-bam")
    new_args["bam_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["bowtie_output"] = Path(new_args["output_prefix"] + "-bowtie")
    new_args["bowtie_output"].mkdir(parents=True, exist_ok=True)
    
    new_args["bowtie_index"] = Path(new_args["bowtie_index"])
    
    new_args["input"] = Path(new_args["input"])
    
    # sequence args
    new_args["irl"] = Seq(new_args["irl"])
    new_args["irr"] = Seq(new_args["irr"])
    new_args["adapter"] = Seq(new_args["adapter"])
    
    # int args
    int_opts = ["verbose", "ntask", "npara"]
    for opts in int_opts:
        new_args[opts] = int(new_args[opts])
        
    return new_args

def preprocess_reads(tpn, adaptor, read_f, read_r, mysample_file, ntask, genome_index_dir, bowtie_output) -> None:
    """Process forward and reverse reads ---> trim adaptors ---> map reads"""
    # for read 1, we need to trim off tpn at 5'
    # for read 2, we need to trim off tpn_rc at 3'
    tpn_rc = tpn.reverse_complement()
    adaptor_rc = adaptor.reverse_complement()

    trim_f = read_f.with_stem("trim-" + read_f.name)
    trim_r = read_r.with_stem("trim-" + read_r.name)

    sam_file = mysample_file.with_suffix(".sam")
    bam_file = mysample_file.with_suffix(".bam")
    bowtie_file = bowtie_output / mysample_file.with_suffix(".txt").name
    
    # https://cutadapt.readthedocs.io/en/stable/guide.html#id4
    # --front or -g, -G is for read 2 (reverse)
    # -g is found by regular 5': Full adapter sequence anywhere, Partial adapter sequence at 5’ end, Full adapter sequence at 5’ end
    # --adapter or -a, -A is for read 2 (reverse)
    # -a is found by regular 3: Full adapter sequence anywhere, Partial adapter sequence at 3’ end, Full adapter sequence at 3’ end
    os.system(
        f"cutadapt -j {ntask} --quiet --discard-untrimmed -g {tpn} -G {adaptor_rc} -a {adaptor} -A {tpn_rc} -o {trim_f} -p {trim_r} {read_f} {read_r}"
    )
    os.system(
        f"bowtie2 -p {ntask} --local -x {genome_index_dir} -q -1 {trim_f} -2 {trim_r} -S {sam_file} > {bowtie_file} 2>&1"
    )
    os.system(f"rm {trim_f}")
    os.system(f"rm {trim_r}")

    os.system(f"samtools sort -@ {ntask} -l 9 -o {bam_file} {sam_file} > /dev/null 2>&1")
    os.system(f"samtools index -@ {ntask} {bam_file}")
    os.system(f"rm {sam_file}")
    sys.stdout.flush()
    sys.stderr.flush()

def preprocess_read_helper(iter_args) -> None:
    """helper function for multiprocessing of preprocess_reads()"""
    row, args = iter_args
    
    data_dir = args["data"]
    bam_output_dir = args["bam_output"]
    bowtie_output_dir = args["bowtie_output"]
    genome_index_dir = args["bowtie_index"]
    ntask = args["ntask"]
    irl_tpn = args["irl"]
    irr_tpn = args["irr"]
    adapter = args["adapter"]
    
    mysample = row[0]
    
    irl_F = data_dir / row[1]
    irl_R = data_dir / row[2]
    irl_file = bam_output_dir / (mysample + "_IRL")
    preprocess_reads(irl_tpn, adapter, irl_F, irl_R, irl_file, ntask, genome_index_dir, bowtie_output_dir)
    
    irr_F = data_dir / row[3]
    irr_R = data_dir / row[4]
    irr_file = bam_output_dir / (mysample + "_IRR")
    preprocess_reads(irr_tpn, adapter, irr_F, irr_R, irr_file, ntask, genome_index_dir, bowtie_output_dir)

def main() -> None:
    main_args = load_args()
    files_df = pd.read_csv(main_args["input"], sep="\t", header=None)
    iter_args = tqdm([ (row[1], main_args) for row in files_df.itertuples() ])
    with Pool(main_args["npara"]) as p:
        [ x for x in p.imap_unordered(preprocess_read_helper, iter_args) ]
        p.close()
    
if __name__ == "__main__":
    main()
