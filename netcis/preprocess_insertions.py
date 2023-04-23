from pathlib import Path
import sys
from multiprocessing import Pool

import pysam
import pandas as pd
import numpy as np


# TODO: changing chromosome names.
# This likely shouldn't stay in here, however, there needs to be someway to state what contigs to keep
# should we expect that contigs are already in a format desireable by the user?
# If that is too advanced for them, then running this tool may be beyond their current skill set

# GRCm39 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/
# https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=&hgsid=1560703641_1YwiSDzyFEZ8nuDrTobTnwtYvReT
chr_dict = {
    "NC_000067.7": "chr1",
    "NC_000068.8": "chr2",
    "NC_000069.7": "chr3",
    "NC_000070.7": "chr4",
    "NC_000071.7": "chr5",
    "NC_000072.7": "chr6",
    "NC_000073.7": "chr7",
    "NC_000074.7": "chr8",
    "NC_000075.7": "chr9",
    "NC_000076.7": "chr10",
    "NC_000077.7": "chr11",
    "NC_000078.7": "chr12",
    "NC_000079.7": "chr13",
    "NC_000080.7": "chr14",
    "NC_000081.7": "chr15",
    "NC_000082.7": "chr16",
    "NC_000083.7": "chr17",
    "NC_000084.7": "chr18",
    "NC_000085.7": "chr19",
    "NC_000086.8": "chrX",
    "NC_000087.8": "chrY",
    "NC_005089.1": "chrM",
}


def get_insertion_properties(insertion, chrdict) -> pd.DataFrame:
    """
    record the insertions stats (direction, +/-, and all that)
    NOTE: here is where additional statistics and or properties for each insertion site can be added
    """
    res = {
        "chr": [chrdict[insertion.reference_name]],
        "pos": [insertion.reference_start],  # 0-based left most coordinate
        "strand +": [insertion.is_forward],
        "ref length": [insertion.reference_length],
        "query length": [
            insertion.infer_query_length()
        ],  # does not include hard-clipped bases
        "read length": [
            insertion.infer_read_length()
        ],  # does include hard-clipped bases. should be equal to len(query_sequence)
        "mapping quality": [insertion.mapping_quality],  # MAPQ: MAPping Quality.
        # MAPQ equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer.
        # A value 255 indicates that the mapping quality is not available.
        # otherwise, the higher the number, the more confident of the quality of the mapping
        # see solution for x in wolfram
        #       254 = -10 * log10(x)
        #       11 = -10 * log10(x)
        "read name": [insertion.query_name],
    }
    res = pd.DataFrame.from_dict(res)
    return res


def read_is_quality(read, is_irr, chr_dict) -> bool:
    # that is paired
    if not read.is_paired:
        return False

    # this is mapped
    if not read.is_mapped:
        return False

    # and has a contig (chromosome) is the predefined dict
    if read.reference_name not in chr_dict.keys():
        return False

    # check if read is forward (+) or reverse (-), then see if 'TA' is present with respects to IRR/IRL orientation
    if read.is_forward:
        if is_irr:  # +, IRR, then starts with TA
            if read.get_forward_sequence()[:2] == "TA":
                return True
        else:  # +, IRL, then ends with TA
            if read.get_forward_sequence()[:-2] == "TA":
                return True
    else:
        if is_irr:  # -, IRR, then ends with TA
            if read.get_forward_sequence()[:-2] == "TA":
                return True
        else:  # -, IRL, then starts with TA
            if read.get_forward_sequence()[:2] == "TA":
                return True

    return False


def process_bam(file, chr_dict, is_irr) -> pd.DataFrame | None:
    """
    Filter out low quality insertions
    This only can run on paired read sequencing data
    is_irr will dictate in what orientation the transposon was sequenced. IRR == forward
    """

    bam = pysam.AlignmentFile(file, "rb")
    insertions = []
    for read1 in bam.fetch():  # multiple_iterators=True
        # only look at read 1
        if not read1.is_read1:
            continue
        # if the read1 is a quality read, then get the insertions properties
        if read_is_quality(read1, is_irr, chr_dict):
            insert_properties = get_insertion_properties(read1, chr_dict)
            insertions.append(insert_properties)
        # check if read 2 (the mate read) is quality and can be used for insertion properties
        else:  # TODO: could be here to check if both reads are quality and match in their reads, then make is a super duper read!
            # must have a mate read that is mapped for .mate() to return properly
            if read1.mate_is_unmapped or (not read1.is_paired):
                continue

            read2 = bam.mate(read1)

            # also check if read2 and read1 mapped to the same reference_name
            if read1.reference_name != read2.reference_name:
                continue

            # then check if the read2 is a quality read and get the insertion properties
            if read_is_quality(read2, is_irr, chr_dict):
                insert_properties = get_insertion_properties(read2, chr_dict)
                insertions.append(insert_properties)

    bam.close()
    # check if there were any inseritons at all to avoid errors from pandas.concat()
    if len(insertions) == 0:
        return None

    insertions_df = pd.concat(insertions, axis=0).reset_index(drop=True)

    # set transposon promoter orientation depending on sequencing library
    # For IRR: + if forward, - if not. For IRL this is reversed
    # using 'strand +' as a cool is easy to change, but it could also have been set this way in get_insertion_properties()
    if is_irr:
        insertions_df["tpn promoter orient +"] = insertions_df["strand +"]
    else:
        insertions_df["tpn promoter orient +"] = ~insertions_df["strand +"]

    # make the orientations easier to read (+/-)
    insertions_df["strand"] = np.where(insertions_df["strand +"], "+", "-")
    insertions_df["tpn promoter orient"] = np.where(
        insertions_df["tpn promoter orient +"], "+", "-"
    )
    insertions_df = insertions_df.drop(["strand +", "tpn promoter orient +"], axis=1)
    return insertions_df

    # # Get the counts by grouping reads that occur at the same chr, pos, strand and tpn promotoer orient
    # group_cols = ['chr', 'pos', 'strand', 'tpn promoter orient']
    # tmp = insertions_df.groupby(by=group_cols, as_index=False, dropna=False)['name'].count()
    # tmp1 = tmp[group_cols]
    # tmp1['count'] = tmp[['name']]
    # # keep track of individual read names to ensure uniqueness of insertion sites
    # tmp1['read names'] = insertions_df.groupby(by=group_cols, dropna=False)['name'].apply(list).reset_index(drop=True)

    # # get the mean of ref, query, read lengths and quality
    # tmp2 = insertions_df.groupby(by=group_cols, as_index=False, dropna=False).mean(numeric_only=True)
    # # change the column names to reflect the mean
    # tmp2 = tmp2.rename({'ref length': 'ref length (mean)', 'query length': 'query length (mean)', 'read length': 'read length (mean)',
    #                   'mapping quality': 'mapping quality (mean)'}, axis=1)

    # # also get median and stdev for mapping quality
    # tmp2['mapping quality (median)'] = insertions_df.groupby(by=group_cols, as_index=False, dropna=False).median(numeric_only=True)['mapping quality']
    # tmp2['mapping quality (stdev)'] = insertions_df.groupby(by=group_cols, as_index=False, dropna=False).std(numeric_only=True)['mapping quality']

    # res_df = tmp1.merge(tmp2, on=group_cols)
    # tmp_names = res_df.pop('read names')
    # res_df.insert(len(res_df.columns.values), 'read names', tmp_names)
    # return res_df


def process_bam_helper(iter_args) -> None:
    row, bam_dir, insertions_dir = iter_args
    mysample = row[0]
    irl_bam = bam_dir / (mysample + "_IRL.bam")
    irr_bam = bam_dir / (mysample + "_IRR.bam")

    # find quality insertion in IRR and IRL libraries and convert them to single insertion site format
    inserts_irl_df = process_bam(file=irl_bam, chr_dict=chr_dict, is_irr=False)
    if (inserts_irl_df is not None):  # if no insertions present, process_bam returns None
        inserts_irl_df["seq library"] = "IRL"

    inserts_irr_df = process_bam(file=irr_bam, chr_dict=chr_dict, is_irr=True)
    if (inserts_irr_df is not None):  # if no insertions present, process_bam returns None
        inserts_irr_df["seq library"] = "IRR"

    # concat of a dataframe and if None just results in the original dataframe
    inserts_df = pd.concat([inserts_irl_df, inserts_irr_df], ignore_index=True)
    # inserts_df = inserts_df.sort_values(['chr', 'pos'], ignore_index=True)

    # # get seq library specific counts
    # count_irr = np.where(inserts_df['seq library'] == 'IRR', inserts_df['count'], 0)
    # count_irl = np.where(inserts_df['seq library'] == 'IRL', inserts_df['count'], 0)
    # inserts_df.insert(6, "count_irr", count_irr)
    # inserts_df.insert(7, "count_irl", count_irl)
    # tmp_read_name = inserts_df.pop('read names')
    # inserts_df.insert(len(inserts_df.columns.values), "read names", tmp_read_name)

    # verify that insertions did not count both read1 and read2
    # do this by checking that the length of 'read names'is the same number as the length of unique read names
    read_names = inserts_df["read name"].to_numpy()
    assert len(np.unique(read_names)) == len(read_names)

    # # show total irr and irl conuts
    # print(f"irr insertions: {inserts_df['count_irr'].sum()}")
    # print(f"irl insertions: {inserts_df['count_irl'].sum()}")

    # save insertions
    inserts_df.to_csv(insertions_dir / (mysample + ".csv"), index=False)


def main() -> None:
    # TODO: change this to load_args using docopt
    output_prefix = sys.argv[1]
    npara = int(sys.argv[2])
    input_files = sys.argv[3]

    bam_dir = Path(output_prefix + "-bam")

    insertions_dir = Path(output_prefix + "-insertions")
    insertions_dir.mkdir(parents=True, exist_ok=True)

    files_df = pd.read_csv(input_files, sep="\t", header=None)
    iter_args = [ (row[1], bam_dir, insertions_dir) for row in files_df.iterrows() ]
    with Pool(npara) as p:
        [ x for x in p.imap_unordered(process_bam_helper, iter_args) ]


if __name__ == "__main__":
    main()
