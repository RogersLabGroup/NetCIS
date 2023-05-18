from pathlib import Path
import sys
import os
from multiprocessing import Pool

import pandas as pd
from Bio.Seq import Seq
from tqdm import tqdm


# TODO: this is hard-coded, can either keep this or give command line parameters for this...
# input sequences are 5' - 3' by standard
irr_tpn = Seq("GGATTAAATGTCAGGAATTGTGAAAA")
irl_tpn = Seq("AAATTTGTGGAGTAGTTGAAAAACGA")
adaptor = Seq("TACCCATACGACGTCCCAGA")


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
    
    # TODO: make sure cutadapt is cutting tags the way we want it to. Write this up in manuscript very specifically
    # https://cutadapt.readthedocs.io/en/stable/guide.html#id4
    # --front or -g, -G is for read 2 (reverse)
    # -g is found by regular 5': Full adapter sequence anywhere, Partial adapter sequence at 5’ end, Full adapter sequence at 5’ end
    # --adapter or -a, -A is for read 2 (reverse)
    # -a is found by regular 3: Full adapter sequence anywhere, Partial adapter sequence at 3’ end, Full adapter sequence at 3’ end
    os.system(
        f"cutadapt -j {ntask} --quiet --discard-untrimmed -g {tpn} -G {adaptor_rc} -a {adaptor} -A {tpn_rc} -o {trim_f} -p {trim_r} {read_f} {read_r}"
        # TODO: look at this other cutadapt again once I have bowtie analysis set
        # cutadapt -a ^FWDPRIMER...RCREVPRIMER -A ^REVPRIMER...RCFWDPRIMER --discard-untrimmed -o out.1.fastq.gz -p out.2.fastq.gz in.1.fastq.gz in.2.fastq.gz
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
    row, data_dir, bam_output_dir, bowtie_output_dir, ntask, genome_index_dir = iter_args
    mysample = row[0]
    irl_F = data_dir / row[1]
    irl_R = data_dir / row[2]
    irr_F = data_dir / row[3]
    irr_R = data_dir / row[4]
    
    irl_file = bam_output_dir / (mysample + "_IRL")
    preprocess_reads(irl_tpn, adaptor, irl_F, irl_R, irl_file, ntask, genome_index_dir, bowtie_output_dir)
    irr_file = bam_output_dir / (mysample + "_IRR")
    preprocess_reads(irr_tpn, adaptor, irr_F, irr_R, irr_file, ntask, genome_index_dir, bowtie_output_dir)

def main() -> None:
    # TODO: change this to load_args using docopt
    data_dir = Path(sys.argv[1])
    output_prefix = sys.argv[2]
    genome_index_dir = Path(sys.argv[3])
    ntask = int(sys.argv[4])  # number of threads to use for bowtie2/samtools
    npara = int(sys.argv[5])  # number of parallel processes for each sample
    input_files = sys.argv[6]

    bam_output_dir = Path(output_prefix + "-bam")
    bam_output_dir.mkdir(parents=True, exist_ok=True)
    
    bowtie_output_dir = Path(output_prefix + "-bowtie")
    bowtie_output_dir.mkdir(parents=True, exist_ok=True)

    files_df = pd.read_csv(input_files, sep="\t", header=None)
    iter_args = [ (row[1], data_dir, bam_output_dir, bowtie_output_dir, ntask, genome_index_dir) for row in files_df.iterrows() ]
    iter_args = tqdm(iter_args)
    with Pool(npara) as p:
        res = [ x for x in p.imap_unordered(preprocess_read_helper, iter_args) ]
    
if __name__ == "__main__":
    main()
