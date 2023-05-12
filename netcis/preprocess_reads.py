from pathlib import Path
import sys
import os
from multiprocessing import Pool

import pandas as pd


# TODO: this is hard-coded, can either keep this or give command line parameters for this...
# I think I might keep it hard-coded, as it is easy for inexperienced users to misspell input args

# TODO: could also only take in the two forward adaptors and use biopython to get the reverse compliment
# primers used to generate libraries
# IRL primer - adaptor1
# IRR primer - adaptor1
# so this is 3 unique tags, but need a totalt of 6 including reverse compliments of tags
# Write this up in manuscript very specificaly
forward_5_adaptor = "GTATGTAAACTTCCGACTTCAACTG"  # adaptor
reverse_5_adaptor = "GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC"  # trnp
forward_3_adaptor = "GTCCCTTAAGCGGAGCCCTATAGTGAGTCGTATTAC"  # trnp
reverse_3_adaptor = "CAGTTGAAGTCGGAAGTTTACATAC"  # adaptor


def preprocess_reads(read_f, read_r, mysample_file, ntask, genome_index_dir) -> None:
    """Process forward and reverse reads ---> trim adaptors ---> map reads"""

    # append temp names to the trimmed files
    trim_f1 = read_f.with_stem("5trim-" + read_f.name)
    trim_r1 = read_r.with_stem("5trim-" + read_r.name)

    trim_f2 = read_f.with_stem("53trim-" + read_f.name)
    trim_r2 = read_r.with_stem("53trim-" + read_r.name)

    trim_f3 = read_f.with_stem("trim-" + read_f.name)
    trim_r3 = read_r.with_stem("trim-" + read_r.name)

    sam_file = mysample_file.with_suffix(".sam")
    bam_file = mysample_file.with_suffix(".bam")

    # TODO: make sure cutadapt is cutting tags the way we want it to. Write this up in manuscript very specifically
    # https://cutadapt.readthedocs.io/en/stable/guide.html#id4
    # --front or -g
    # -G is for read 2 (reverse)
    # # -g is found by regular 5': Full adapter sequence anywhere, Partial adapter sequence at 5’ end, Full adapter sequence at 5’ end
    # # -g ^ is found by anchored 5': Full adapter sequence at 5’ end
    os.system(
        f"cutadapt --quiet -j {ntask} --discard-untrimmed -g {forward_5_adaptor} -o {trim_f1} -p {trim_r1} {read_f} {read_r}"
    )
    os.system(
        f"cutadapt --quiet -j {ntask} -G ^{reverse_5_adaptor} -o {trim_f2} -p {trim_r2} {trim_f1} {trim_r1}"
    )
    os.system(f"rm {trim_f1}")
    os.system(f"rm {trim_r1}")

    # --adapter or -a
    # -A is for read 2 (reverse)
    # -a is found by regular 3: Full adapter sequence anywhere, Partial adapter sequence at 3’ end, Full adapter sequence at 3’ end
    os.system(
        f"cutadapt --quiet -j {ntask} -a {forward_3_adaptor} -A {reverse_3_adaptor} -o {trim_f3} -p {trim_r3} {trim_f2} {trim_r2}"
    )
    os.system(f"rm {trim_f2}")
    os.system(f"rm {trim_r2}")

    os.system(
        f"bowtie2 -p {ntask} --local --quiet -x {genome_index_dir} -q -1 {trim_f3} -2 {trim_r3} -S {sam_file}"
    )
    os.system(f"rm {trim_f3}")
    os.system(f"rm {trim_r3}")

    os.system(f"samtools sort -@ {ntask} -m 4G -l 9 -o {bam_file} {sam_file}")
    os.system(f"samtools index -@ {ntask} {bam_file}")
    os.system(f"rm {sam_file}")


def preprocess_read_helper(iter_args) -> None:
    row, data_dir, bam_output_dir, ntask, genome_index_dir = iter_args
    mysample = row[0]
    irl_F = data_dir / row[1]
    irl_R = data_dir / row[2]
    irr_F = data_dir / row[3]
    irr_R = data_dir / row[4]

    irl_file = bam_output_dir / (mysample + "_IRL")
    preprocess_reads(irl_F, irl_R, irl_file, ntask, genome_index_dir)

    irr_file = bam_output_dir / (mysample + "_IRR")
    preprocess_reads(irr_F, irr_R, irr_file, ntask, genome_index_dir)


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

    files_df = pd.read_csv(input_files, sep="\t", header=None)
    iter_args = [ (row[1], data_dir, bam_output_dir, ntask, genome_index_dir) for row in files_df.iterrows() ]
    with Pool(npara) as p:
        [ x for x in p.imap_unordered(preprocess_read_helper, iter_args) ]


if __name__ == "__main__":
    main()
