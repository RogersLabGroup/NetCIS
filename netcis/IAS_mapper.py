#!/usr/bin/python
from pathlib import Path
import sys
import os
import itertools


def file_len(fname):
    with open(fname, 'r') as f:
        return len(f.readlines())

def cigar_caller(mystring):
    mynumb = ""
    total = 0
    for char in mystring:
        if char.isdigit() == True:
            mynumb += char
        elif char == "M":
            total += int(mynumb)
            mynumb = ""
        elif char == "I":
            total -= int(mynumb)
            mynumb = ""
        elif char == "D":
            total += int(mynumb)*2
            mynumb = ""
    return total
    
def add_results(result, tag, ind_1, ind_2, max_count):
    mapcheck = False
    if tag in result:
        result[tag][ind_1] += 1  # FIXME: 
        if result[tag][ind_1] > max_count:  # FIXME: 
            max_count = result[tag][ind_1]  # FIXME: 
        mapcheck = True
    else:
        result[tag] = [ind_2, ind_1]  # FIXME: 
        mapcheck = True
    return result, max_count, mapcheck

def add_ligation(ligation, tag, ind_1, length_read):
    if tag in ligation:
        if length_read not in ligation[tag][ind_1]:  # FIXME: 
            ligation[tag][ind_1].append(length_read)  # FIXME: 
    else:
        ligation[tag] = [[], []]
        ligation[tag][ind_1].append(length_read)  # FIXME:
    return ligation
       
def samparse(sam_file, read_type, results, ligations):
    # read_type can be IRL or IRR and denotes the type of sam file we are working on
    if read_type == 'IRL':
        ind1 = 0
        ind2 = 1
        fwd_strand = '-'
        rev_strand = '+'
    elif read_type == 'IRR':
        ind1 = 1
        ind2 = 0
        fwd_strand = '+'
        rev_strand = '-'
    else:
        return
    
    read1_flags = [69, 73, 89, 99, 81, 83, 89]
    read1_flags_pos = [69, 73, 89, 99]
    read1_flags_neg = [81, 83, 89]
    
    read2_flags = [161, 163, 169, 145, 147]
    read2_flags_pos = [161, 163, 169]
    read2_flags_neg = [145, 147]
    
    max_read = 0
    with open(sam_file, "r") as SAM:
        count = 0
        linecount = file_len(sam_file)
        for line, line2 in itertools.zip_longest(SAM, SAM, fillvalue=''):
            if linecount - count > 2:
                count += 2
                mylist = line.split("\t")
                mylist2 = line2.split("\t")
                flag1 = int(mylist[1])
                flag2 = int(mylist2[1])
                
                if flag1 in read1_flags:
                    myflag = flag1
                    chrom = mylist[2]
                    address = int(mylist[3])
                    length = cigar_caller(mylist[5])
                    myseq = mylist[9]
                elif flag2 in read1_flags:
                    myflag = flag2
                    chrom = mylist2[2]
                    address = int(mylist2[3])
                    length = cigar_caller(mylist2[5])
                    myseq = mylist2[9]
                else:
                    length = 0  # or just continue to the next pair of reads
                    continue
                    
                    
                # Runs if first line contains the mapping data for the transposon read
                # Executes if transposon IRL read mapped to + strand
                if length >= 12:
                    if myflag in read1_flags_pos:
                        if myseq[0:2] == "TA":
                            mytag = f"{chrom}:{address}{fwd_strand}"  # FIXME: flips strand for IRL (-) or for IRR (+)
                            results, max_read, mapcheck = add_results(results, mytag, ind1, ind2, max_read)
                            ligations = add_ligation(ligations, mytag, ind1, length) 
                    elif myflag in read1_flags_neg:
                        if myseq[-2:] == "TA":
                            mytag = f"{chrom}:{address+(length-2)}{rev_strand}" # FIXME: flips strand for IRL (+) or for IRR (-)
                            results, max_read, mapcheck = add_results(results, mytag, ind1, ind2, max_read)
                            ligations = add_ligation(ligations, mytag, ind1, length)
                        
                    if mapcheck == False:  # paired read only if the transposon read cannot be mapped precisely
                        if flag1 in read2_flags:
                            myflag = flag1
                            chrom =  mylist[2]
                            address = int(mylist[3])
                            length = cigar_caller(mylist[5])
                            myseq = mylist[9]
                        elif flag2 in read2_flags:
                            myflag = flag2
                            chrom = mylist2[2]
                            address = int(mylist2[3])
                            length = cigar_caller(mylist2[5])
                            myseq = mylist2[9]
                            
                        if length >= 10:
                            if myflag in read2_flags_pos:
                                if myseq[0:2] == "TA":  # this line is the same for either IRL or IRR file
                                    mytag = f"{chrom}:{address}-"  # this line is the same for either IRL or IRR file
                                    results, max_read, mapcheck = add_results(results, mytag, ind1, ind2, max_read)
                                    ligations = add_ligation(ligations, mytag, ind1, length)
                                    
                            elif myflag in read2_flags_neg:
                                if myseq[-2:] == "TA": # this line is the same for either IRL or IRR file
                                    mytag = f"{chrom}:{address+(length-2)}+"  # flips strand since it's the IRL read
                                    results, max_read, mapcheck = add_results(results, mytag, ind1, ind2, max_read)
                                    ligations = add_ligation(ligations, mytag, ind1, length)            
    return results, ligations, max_read


def samparse_orig(irlsam, irrsam, output_dir):
    results = {}
    ligation = {}
    read1_flags = [69,73,89,99,81,83,89]
    read2_flags = [161,163,169,145,147]
    read1_flags_pos = [69,73,89,99]
    read1_flags_neg = [81,83,89]
    read2_flags_pos = [161,163,169]
    read2_flags_neg = [145,147]
    IRLmax,IRRmax = 0,0
    
    with open(irlsam, "r") as SAM:
        count = 0
        linecount = file_len(irlsam)
        for line, line2 in itertools.zip_longest(SAM,SAM,fillvalue=''):
            if linecount - count > 2:
                count += 2
                mapcheck = False
                mylist = line.split("\t")
                mylist2 = line2.split("\t")
                flag1,flag2 = int(mylist[1]), int(mylist2[1])
                if flag1 in read1_flags:
                    myflag = flag1
                    chrom, address = mylist[2], int(mylist[3])
                    length = cigar_caller(mylist[5])
                    myseq = mylist[9]
                elif flag2 in read1_flags:
                    myflag = flag2
                    chrom, address = mylist2[2], int(mylist2[3])
                    length = cigar_caller(mylist2[5])
                    myseq = mylist2[9]
                else:
                    length = 0
                #Runs if first line contains the mapping data for the transposon read
                #Executes if transposon IRL read mapped to + strand
                if length >= 12:  # there is a minimum len input for hisat2, so it has to be larger than that
                    if myflag in read1_flags_pos:
                        if myseq[0:2] == "TA":
                            mytag = chrom +":"+ str(address) +"-" #flips strand for IRL
                            if mytag in results:
                                results[mytag][0] += 1
                                if results[mytag][0] > IRLmax:
                                    IRLmax = results[mytag][0]
                                mapcheck = True
                            else:
                                results[mytag] = [1,0]
                                mapcheck = True
                            if mytag in ligation:
                                if length not in ligation[mytag][0]:
                                    ligation[mytag][0].append(length)
                            else:
                                ligation[mytag] = [[],[]]
                                ligation[mytag][0].append(length)
                    elif myflag in read1_flags_neg:
                        if myseq[-2:] == "TA":
                            mytag = chrom +":"+ str(address+(length-2)) +"+" #flips strand for IRL
                            if mytag in results:
                                results[mytag][0] += 1
                                if results[mytag][0] > IRLmax:
                                    IRLmax = results[mytag][0]
                                mapcheck = True
                            else:
                                results[mytag] = [1,0]
                                mapcheck = True
                            if mytag in ligation:
                                if length not in ligation[mytag][0]:
                                    ligation[mytag][0].append(length)
                            else:
                                ligation[mytag] = [[],[]]
                                ligation[mytag][0].append(length)
                    if mapcheck == False: #paired read only if the tnp read cannot be mapped precisely
                        if flag1 in read2_flags:
                            myflag = flag1
                            chrom, address = mylist[2], int(mylist[3])
                            length = cigar_caller(mylist[5])
                            myseq = mylist[9]
                        elif flag2 in read2_flags:
                            myflag = flag2
                            chrom, address = mylist2[2], int(mylist2[3])
                            length = cigar_caller(mylist2[5])
                            myseq = mylist2[9]
                        if length >= 10:
                            if myflag in read2_flags_pos:
                                if myseq[0:2] == "TA":
                                    mytag = chrom +":"+ str(address) +"-"
                                    if mytag in results:
                                        results[mytag][0] += 1
                                        if results[mytag][0] > IRLmax:
                                            IRLmax = results[mytag][0]
                                        mapcheck = True
                                    else:
                                        results[mytag] = [1,0]
                                        mapcheck = True
                                    if mytag in ligation:
                                        if length not in ligation[mytag][0]:
                                            ligation[mytag][0].append(length)
                                    else:
                                        ligation[mytag] = [[],[]]
                                        ligation[mytag][0].append(length)
                            elif myflag in read2_flags_neg:
                                if myseq[-2:] == "TA":
                                    mytag = chrom +":"+ str(address+(length-2)) +"+" #flips strand since it's the IRL read
                                    if mytag in results:
                                        results[mytag][0] += 1
                                        if results[mytag][0] > IRLmax:
                                            IRLmax = results[mytag][0]
                                        mapcheck = True
                                    else:
                                        results[mytag] = [1,0]
                                        mapcheck = True
                                    if mytag in ligation:
                                        if length not in ligation[mytag][0]:
                                            ligation[mytag][0].append(length)
                                    else:
                                        ligation[mytag] = [[],[]]
                                        ligation[mytag][0].append(length)
    with open(irrsam, "r") as SAM:
        count = 0
        linecount = file_len(irrsam)
        for line, line2 in itertools.zip_longest(SAM,SAM,fillvalue=''):
            if linecount - count > 2:
                count += 2
                mapcheck = False
                mylist = line.split("\t")
                mylist2 = line2.split("\t")
                flag1, flag2 = int(mylist[1]), int(mylist2[1])
                if flag1 in read1_flags:
                    myflag = flag1
                    chrom, address = mylist[2], int(mylist[3])
                    length = cigar_caller(mylist[5])
                    myseq = mylist[9]
                elif flag2 in read1_flags:
                    myflag = flag2
                    chrom, address = mylist2[2], int(mylist2[3])
                    length = cigar_caller(mylist2[5])
                    myseq = mylist2[9]
                else:
                    length = 0
                #Runs if first line contains the mapping data for the transposon read
                #Executes if transposon IRL read mapped to + strand
                if length >= 12:
                    if myflag in read1_flags_pos:
                        if myseq[0:2] == "TA":
                            mytag = chrom +":"+ str(address) +"+" #same strand for IRR
                            if mytag in results:
                                results[mytag][1] += 1
                                if results[mytag][1] > IRRmax:
                                    IRRmax = results[mytag][1]
                                mapcheck = True
                            else:
                                results[mytag] = [0,1]
                                mapcheck = True
                            if mytag in ligation:
                                if length not in ligation[mytag][1]:
                                    ligation[mytag][1].append(length)
                            else:
                                ligation[mytag] = [[],[]]
                                ligation[mytag][1].append(length)
                    elif myflag in read1_flags_neg:
                        if myseq[-2:] == "TA":
                            mytag = chrom +":"+ str(address+(length-2)) +"-" #flips strand for IRR
                            if mytag in results:
                                results[mytag][1] += 1
                                if results[mytag][1] > IRRmax:
                                    IRRmax = results[mytag][1]
                                mapcheck = True
                            else:
                                results[mytag] = [0,1]
                                mapcheck = True
                            if mytag in ligation:
                                if length not in ligation[mytag][1]:
                                    ligation[mytag][1].append(length)
                            else:
                                ligation[mytag] = [[],[]]
                                ligation[mytag][1].append(length)
                    if mapcheck == False: #paired read only analyzed if the transposon read cannot be mapped precisely
                        if flag1 in read2_flags:
                            myflag = flag1
                            chrom, address = mylist[2], int(mylist[3])
                            length = cigar_caller(mylist[5])
                            myseq = mylist[9]
                        elif flag2 in read2_flags:
                            myflag = flag2
                            chrom, address = mylist2[2], int(mylist2[3])
                            length = cigar_caller(mylist2[5])
                            myseq = mylist2[9]
                        if length >= 10:
                            if myflag in read2_flags_pos:
                                if myseq[0:2] == "TA":
                                    mytag = chrom +":"+ str(address) +"-"
                                    if mytag in results:
                                        results[mytag][1] += 1
                                        if results[mytag][1] > IRRmax:
                                            IRRmax = results[mytag][1]
                                        mapcheck = True
                                    else:
                                        results[mytag] = [0,1]
                                        mapcheck = True
                                    if mytag in ligation:
                                        if length not in ligation[mytag][1]:
                                            ligation[mytag][1].append(length)
                                    else:
                                        ligation[mytag] = [[],[]]
                                        ligation[mytag][1].append(length)
                            elif myflag in read2_flags_neg:
                                if myseq[-2:] == "TA":
                                    mytag = chrom +":"+ str(address+(length-2)) +"+" #flips strand since it's the IRL read
                                    if mytag in results:
                                        results[mytag][1] += 1
                                        if results[mytag][1] > IRRmax:
                                            IRRmax = results[mytag][1]
                                        mapcheck = True
                                    else:
                                        results[mytag] = [0,1]
                                        mapcheck = True
                                    if mytag in ligation:
                                        if length not in ligation[mytag][1]:
                                            ligation[mytag][1].append(length)
                                    else:
                                        ligation[mytag] = [[],[]]
                                        ligation[mytag][1].append(length)
    
    myfile = irlsam.name[:-8]
    output = open(output_dir / (myfile + ".uniq"), "w")
    print(len(results.keys()))
    for key in results:
        mylist = key.split(":")
        chrom,address = mylist[0],mylist[1][:-1]
        ligationlist = ligation[key]
        IRLlp = len(ligationlist[0])
        IRRlp = len(ligationlist[1])
        strand = mylist[1][-1:]
        myvals = results[key]
        IRL,IRR = myvals[0],myvals[1]
        try:
            IRLnorm = round((float(IRL)/IRLmax*100),3)  
            # clonality of the tumor cells, or rather the ordering of the insertion events of the tumor 
            # TODO: also need to quanitfy the read quality of a tumor, overall low quality reads may need more analysis to determine actual quality
        except:
            IRLnorm = 0
        try:
            IRRnorm = round((float(IRR)/IRRmax*100),3)
        except:
            IRRnorm = 0
        mystring = myfile +"\t"+ chrom +"\t"+ str(address) +"\t"+  str(IRL) +"\t"+ str(IRR) +"\t"+ strand +"\t"+ str(IRLlp) +"\t"+ str(IRRlp) +"\t"+ str(IRLnorm) +"\t"+ str(IRRnorm) +"\n"
        output.write(mystring)
    output.close()
    
    
def main():
    data_dir = Path(sys.argv[1])
    sam_output_dir = Path(sys.argv[2])
    uniq_output_dir = Path(sys.argv[3])
    genome_index_dir = Path(sys.argv[4])
    N = int(sys.argv[5])  # number of threds to use
    
    mysample = sys.argv[6]
    irl_F = data_dir / Path(sys.argv[7])
    irl_R = data_dir / Path(sys.argv[8])
    irr_F = data_dir / Path(sys.argv[9])
    irr_R = data_dir / Path(sys.argv[10])
    

    # Process IRL reads ---> trim adaptor ---> map reads
    # append temp names to the files
    trim_irl_f1 = data_dir / (irl_F.name.rstrip("".join(irl_F.suffixes)) + "-5trim" + "".join(irl_F.suffixes))
    trim_irl_r1 = data_dir / (irl_R.name.rstrip("".join(irl_R.suffixes)) + "-5trim" + "".join(irl_R.suffixes))
    trim_irl_f2 = data_dir / (irl_F.name.rstrip("".join(irl_F.suffixes)) + "-53trim" + "".join(irl_F.suffixes))
    trim_irl_r2 = data_dir / (irl_R.name.rstrip("".join(irl_R.suffixes)) + "-53trim" + "".join(irl_R.suffixes))
    trim_irl_f3 = data_dir / (irl_F.name.rstrip("".join(irl_F.suffixes)) + "-trim" + "".join(irl_F.suffixes))
    trim_irl_r3 = data_dir / (irl_R.name.rstrip("".join(irl_R.suffixes)) + "-trim" + "".join(irl_R.suffixes))
    irl_sam = sam_output_dir / (mysample + "_IRL.sam")
    
    # TODO: change cutadapt to latest version
    os.system(f"~/.local/bin/cutadapt --quiet -j {N} --discard-untrimmed -g GTATGTAAACTTCCGACTTCAACTG -o {trim_irl_f1} -p {trim_irl_r1} {irl_F} {irl_R}")

    os.system(f"~/.local/bin/cutadapt --quiet -j {N} -G ^GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC -o {trim_irl_f2} -p {trim_irl_r2} {trim_irl_f1} {trim_irl_r1}")
    os.system(f"rm {trim_irl_f1}")
    os.system(f"rm {trim_irl_r1}")
    
    os.system(f"~/.local/bin/cutadapt --quiet -j {N} -a GTCCCTTAAGCGGAGCCCTATAGTGAGTCGTATTAC -A CAGTTGAAGTCGGAAGTTTACATAC -o {trim_irl_f3} -p {trim_irl_r3} {trim_irl_f2} {trim_irl_r2}")
    os.system(f"rm {trim_irl_f2}")
    os.system(f"rm {trim_irl_r2}")

    # python edited_IAS_mapper.py example_fastq/IAS_input-hisat2.tsv 
    # os.system(f"hisat2 -p {N} --no-hd --no-sq -x {genome_index_dir} -q -1 {trim_irl_f3} -2 {trim_irl_r3} -S {irl_sam}")  # 0.43% overall alignment rate
    
    # python edited_IAS_mapper.py example_fastq/IAS_input-bowtie2.tsv
    # os.system(f"bowtie2 -p {N} --no-hd --no-sq -x {genome_index_dir} -q -1 {trim_irl_f3} -2 {trim_irl_r3} -S {irl_sam}")  # 2.77% overall alignment rate
    os.system(f"bowtie2 -p {N} --local --quiet -x {genome_index_dir} -q -1 {trim_irl_f3} -2 {trim_irl_r3} -S {irl_sam}")  # 16.16% overall alignment rate
    os.system(f"samtools view -@8 {irl_sam} -h -b -o {irl_sam.with_suffix('.bam')}")
    os.system(f"rm {trim_irl_f3}")
    os.system(f"rm {trim_irl_r3}")

    
    # hisat2 takes much longer than bowtie2        



    # Process IRR reads ---> trim adaptor ---> map reads
    trim_irr_f1 = data_dir / (irr_F.name.rstrip("".join(irr_F.suffixes)) + "-5trim" + "".join(irr_F.suffixes))
    trim_irr_r1 = data_dir / (irr_R.name.rstrip("".join(irr_R.suffixes)) + "-5trim" + "".join(irr_R.suffixes))
    trim_irr_f2 = data_dir / (irr_F.name.rstrip("".join(irr_F.suffixes)) + "-53trim" + "".join(irr_F.suffixes))
    trim_irr_r2 = data_dir / (irr_R.name.rstrip("".join(irr_R.suffixes)) + "-53trim" + "".join(irr_R.suffixes))
    trim_irr_f3 = data_dir / (irr_F.name.rstrip("".join(irr_F.suffixes)) + "-trim" + "".join(irr_F.suffixes))
    trim_irr_r3 = data_dir / (irr_R.name.rstrip("".join(irr_R.suffixes)) + "-trim" + "".join(irr_R.suffixes))
    irr_sam = sam_output_dir / (mysample + "_IRR.sam")
    
    os.system(f"~/.local/bin/cutadapt --quiet -j {N} --discard-untrimmed -g GTATGTAAACTTCCGACTTCAACTG -o {trim_irr_f1} -p {trim_irr_r1} {irr_F} {irr_R}")
    
    os.system(f"~/.local/bin/cutadapt --quiet -j {N} -G ^GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC -o {trim_irr_f2} -p {trim_irr_r2} {trim_irr_f1} {trim_irr_r1}")
    os.system(f"rm {trim_irr_f1}")
    os.system(f"rm {trim_irr_r1}")
    
    os.system(f"~/.local/bin/cutadapt --quiet -j {N} -a GTCCCTTAAGCGGAGCCCTATAGTGAGTCGTATTAC -A CAGTTGAAGTCGGAAGTTTACATAC -o {trim_irr_f3} -p {trim_irr_r3} {trim_irr_f2} {trim_irr_r2}")
    os.system(f"rm {trim_irr_f2}")
    os.system(f"rm {trim_irr_r2}")
    
    # os.system(f"hisat2 -p {N} --no-hd --no-sq -x {genome_index_dir} -q -1 {trim_irr_f3} -2 {trim_irr_r3} -S {irr_sam}")  # 4.02% overall alignment rate
    # os.system(f"bowtie2 -p {N} --no-hd --no-sq -x {genome_index_dir} -q -1 {trim_irr_f3} -2 {trim_irr_r3} -S {irr_sam}")  # 7.32% overall alignment rate
    os.system(f"bowtie2 -p {N} --local --quiet -x {genome_index_dir} -q -1 {trim_irr_f3} -2 {trim_irr_r3} -S {irr_sam}")  # 40.91% overall alignment rate
    os.system(f"samtools view -@8 {irr_sam} -h -b -o {irr_sam.with_suffix('.bam')}")
    os.system(f"rm {trim_irr_f3}")
    os.system(f"rm {trim_irr_r3}")
    
    # TODO: try minimap2? or something else?
    # create an index first and then map
    # ./minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa
    # ./minimap2 -a MT-human-ont.mmi test/MT-orang.fa > test.sam
    # os.system(f"bowtie2 -p {N} --no-hd --no-sq --local --quiet -x {genome_index_dir} -q -1 {trim_irr_f3} -2 {trim_irr_r3} -S {irr_sam}")  # 40.91% overall alignment rate

    os.system(f"samtools view {irl_sam.with_suffix('.bam')} --no-header -o {irl_sam}")
    os.system(f"samtools view {irr_sam.with_suffix('.bam')} --no-header -o {irr_sam}")
    samparse_orig(irl_sam, irr_sam, uniq_output_dir)
    os.system(f"rm {irl_sam}")
    os.system(f"rm {irr_sam}")
    
    
    
    # # Combine mapped IRL and IRR reads ---> generate .UNIQ file
    # # TODO: there has to be a better way to do this...
    # # TODO: maybe combine all these .uniq files into one tabular file with pandas?
    # results_dict = {}
    # ligations_dict = {}
    # results_dict, ligations_dict, IRL_max = samparse(irl_sam, 'IRL', results_dict, ligations_dict)
    # results_dict, ligations_dict, IRR_max = samparse(irr_sam, 'IRR', results_dict, ligations_dict)
    
    # myfile = irl_sam[:-8]
    # with open(myfile + ".uniq", "w") as output:  # TODO: do we need to append or write to file?
    #     for key in results_dict:
    #         mylist = key.split(":")
    #         chrom,address = mylist[0],mylist[1][:-1]
    #         ligationlist = ligations_dict[key]
    #         IRLlp = len(ligationlist[0])
    #         IRRlp = len(ligationlist[1])
    #         strand = mylist[1][-1:]
    #         myvals = results_dict[key]
    #         IRL = myvals[0]
    #         IRR = myvals[1]
    #         try:
    #             IRLnorm = round((float(IRL) / IRL_max * 100), 3)
    #         except:
    #             IRLnorm = 0
    #         try:
    #             IRRnorm = round((float(IRR) / IRR_max  *100), 3)
    #         except:
    #             IRRnorm = 0
                
    #         # mystring = myfile +"\t"+ chrom +"\t"+ str(address) +"\t"+  str(IRL) +"\t"+ str(IRR) +"\t"+ strand +"\t"+ str(IRLlp) +"\t"+ str(IRRlp) +"\t"+ str(IRLnorm) +"\t"+ str(IRRnorm) +"\n"
    #         # print(mystring, file=output, end="")
    #         print("\t".join([myfile, chrom, address, IRL, IRR, strand, IRLlp, IRRlp, IRLnorm, IRRnorm]), file=output)
    
    # os.system("rm " + irl_sam)
    # os.system("rm " + irr_sam) 
    

if __name__  == "__main__":
    main()
    