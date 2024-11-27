import sys, subprocess, time
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
     -v, --verbose=N                    print more verbose information if available using 0, 1 or 2 [default: 0]
     -t, --ntask=INT                    number of threads to use that will speed up cutadapt, bowtie2, and samtools [default: 1]
     -n, --npara=INT                    number of parallel processes to process multiple files at the same time [default: 1]
     -q, --mapq=INT                     minimum mapping quality score to keep reads [default: 13]   
    """

    # remove "--" from args
    args = { key.split("-")[-1]: value for key, value in docopt(doc).items() }
    
    # int args
    int_opts = ["verbose", "ntask", "npara"]
    for opts in int_opts:
        args[opts] = int(args[opts])
    
    if args["verbose"] > 1:
        print("Arguements given")
        for key, item in args.items():
            print(f"\t{key}: {item}")
        print("\n")
        
    # files and directory args
    args["data"] = Path(args["data"])
    
    args["bam_output"] = Path(args["output_prefix"] + "-bam")
    args["bam_output"].mkdir(parents=True, exist_ok=True)
    
    args["report_output"] = args["bam_output"].with_name("reports-preprocessing")
    args["report_output"].mkdir(parents=True, exist_ok=True)
    
    args["bowtie_index"] = Path(args["bowtie_index"])
    
    args["input"] = Path(args["input"])
    
    # sequence args
    args["irl"] = Seq(args["irl"])
    args["irr"] = Seq(args["irr"])
    args["primer"] = Seq(args["primer"])
    
    return args

def run_command(command):
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"Error in command: {command}")
        print(result.stderr.decode())
    return result

def run_preprocessing(cutadapt, ntask, genome_index_dir, trim1_f, trim1_r, pre_bam_file, bowtie_report, bam_file, bam_report, report_output):
    # trim reads
    res_cutadapt = run_command(cutadapt)
    
    # map reads
    mapper = f"bowtie2 -p {ntask} --very-sensitive-local --local -x {genome_index_dir} -1 {trim1_f} -2 {trim1_r}"
    mapper2 = f"samtools view -h -@ {ntask} -b - > {pre_bam_file} 2> {report_output / bowtie_report}"
    res_bowtie = run_command(f"{mapper} | {mapper2}")
    
    # filter reads
    sam_sort = f"samtools sort -@ {ntask} -o {pre_bam_file} {pre_bam_file}"
    # sam_filter = f"samtools view -h -@ {ntask} -f 3 -1 -o {bam_file} {pre_bam_file}"
    sam_filter = f"samtools view -h -@ {ntask} -q 13 -f 67 -F 12 -1 -o {bam_file} {pre_bam_file}"
    sam_index = f"samtools index -@ {ntask} {bam_file}"
    sam_index2 = f"samtools index -@ {ntask} {pre_bam_file}"
    sam_stats = f"samtools idxstats {bam_file} > {report_output / bam_report}"
    rm_files = f"rm {trim1_f} {trim1_r}"
    res_idxstats = run_command(f"{sam_sort} ; {sam_filter} ; {sam_index} ; {sam_index2} ; {sam_stats} ; {rm_files}")
                   
    sys.stdout.flush()
    sys.stderr.flush()
    
    assert (res_cutadapt.returncode == 0), f"Error in cutadapt\n\n{res_cutadapt.stderr}"
    assert (res_bowtie.returncode == 0), f"Error in bowtie\n\n{res_bowtie.stderr}"
    assert (res_idxstats.returncode == 0), f"Error in idxstats\n\n{res_idxstats.stderr}"

def preprocess_reads(tpn: Seq, primer: Seq, read_f: Path, read_r: Path, mysample_file: Path, library, tpn_orient, ntask, genome_index_dir, report_output):
    """Process forward and reverse reads: trim transposon and primer, map reads, save to bam file"""

    if tpn_orient == "+":
        trim1_f = mysample_file.with_name("trim1-orient_pos-" + read_f.name)
        trim1_r = mysample_file.with_name("trim1-orient_pos-" + read_r.name)
        cutadapt_report = mysample_file.with_suffix(".cutadapt-orient_pos.txt").name
        
        pre_bam_file = mysample_file.with_suffix(".prefiltering-orient_pos.bam")
        bowtie_report = mysample_file.with_suffix(".bowtie-orient_pos.txt").name
        
        bam_file = mysample_file.with_suffix(".orient_pos.bam")
        bam_report = mysample_file.with_suffix(".idxstats-orient_pos.txt").name
        
        if library == "IRL":
            # IRL and tpn in forward orientation with strand
            cutadapt = (
                f"cutadapt -j {ntask} -m 20 --discard-untrimmed --pair-filter=any "
                f"-g ^{primer} -a {tpn} -G ^{tpn.reverse_complement()} -A {primer.reverse_complement()} "
                f"-o {trim1_f} -p {trim1_r} {read_f} {read_r} > {report_output / cutadapt_report}"
            )
        else:
            # IRR and tpn in forward orientation with strand
            cutadapt = (
                f"cutadapt -j {ntask} -m 20 --discard-untrimmed --pair-filter=any "
                f"-g ^{tpn} -a {primer} -G ^{primer.reverse_complement()} -A {tpn.reverse_complement()} "
                f"-o {trim1_f} -p {trim1_r} {read_f} {read_r} > {report_output / cutadapt_report}"
            )
            
    else:
        trim1_f = mysample_file.with_name("trim1-orient_neg-" + read_f.name)
        trim1_r = mysample_file.with_name("trim1-orient_neg-" + read_r.name)
        cutadapt_report = mysample_file.with_suffix(".cutadapt-orient_neg.txt").name
        
        pre_bam_file = mysample_file.with_suffix(".prefiltering-orient_neg.bam")
        bowtie_report = mysample_file.with_suffix(".bowtie-orient_neg.txt").name
        
        bam_file = mysample_file.with_suffix(".orient_neg.bam")
        bam_report = mysample_file.with_suffix(".idxstats-orient_neg.txt").name

        if library == "IRL":
            # IRL and tpn in reverse orientation against strand
            cutadapt = (
                f"cutadapt -j {ntask} -m 20 --discard-untrimmed --pair-filter=any "
                f"-g ^{tpn.reverse_complement()} -a {primer.reverse_complement()} -G ^{primer} -A {tpn} "
                f"-o {trim1_f} -p {trim1_r} {read_f} {read_r} > {report_output / cutadapt_report}"
            )
        else:
            # IRR and tpn in reverse orientation against strand
            cutadapt = (
                f"cutadapt -j {ntask} -m 20 --discard-untrimmed --pair-filter=any "
                f"-g ^{primer.reverse_complement()} -a {tpn.reverse_complement()} -G ^{tpn} -A {primer} "
                f"-o {trim1_f} -p {trim1_r} {read_f} {read_r} > {report_output / cutadapt_report}"
            )
    
    run_preprocessing(
        cutadapt,
        ntask, genome_index_dir, trim1_f, trim1_r, pre_bam_file, bowtie_report, 
        bam_file, bam_report,
        report_output,
        )

def preprocess_read_helper(iter_args) -> None:
    """helper function for multiprocessing of preprocess_reads()"""
    row, library, tpn_orient, args = iter_args
    
    data_dir = args["data"]
    bam_output_dir = args["bam_output"]
    genome_index_dir = args["bowtie_index"]
    ntask = args["ntask"]
    irl_tpn = args["irl"]
    irr_tpn = args["irr"]
    primer = args["primer"]
    report_output_dir = args["report_output"]
    
    # separate IRL and IRR to run in parallel
    mysample = row.iloc[0]
    
    if library == 'IRL':
        irl_F = data_dir / row.iloc[1]
        irl_R = data_dir / row.iloc[2]
        irl_file = bam_output_dir / (mysample + f"_{library}")
        preprocess_reads(irl_tpn, primer, irl_F, irl_R, irl_file, library, tpn_orient, ntask, genome_index_dir, report_output_dir)
    
    else:  # library == 'IRR'
        irr_F = data_dir / row.iloc[3]
        irr_R = data_dir / row.iloc[4]
        irr_file = bam_output_dir / (mysample + f"_{library}")
        preprocess_reads(irr_tpn, primer, irr_F, irr_R, irr_file, library, tpn_orient, ntask, genome_index_dir, report_output_dir)
        
def main() -> None:
    main_args = load_args()
    files_df = pd.read_csv(main_args["input"], sep="\t")
    iter_args = []
    for row in files_df.iterrows():
        # input file row, library, tpn orientation, args
        iter_args.append((row[1], 'IRL', '+', main_args))
        iter_args.append((row[1], 'IRL', '-', main_args))
        iter_args.append((row[1], 'IRR', '+', main_args))
        iter_args.append((row[1], 'IRR', '-', main_args))
        
    if main_args['verbose']:
        print("preprocess_reads.py:")
        iter_args = tqdm(iter_args)

    with Pool(main_args["npara"]) as p:
        [ x for x in p.imap_unordered(preprocess_read_helper, iter_args, chunksize=1) ]
        p.close()
    
if __name__ == "__main__":
    main()
