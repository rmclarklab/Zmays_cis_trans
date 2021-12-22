#!/home/clarklab/Local/anaconda3/bin/python3.7

"""This script is to restore snp and indel in vcf based on reference genome.
SNP variants would show with capitalized characters.
INDEL marked '?' as insertion, and '_' as deletion. 
NOTE: sites with multiple alternatives are neglected. 
From Meiyuan Ji, on 10.29.2019."""

from Bio import SeqIO
import re
import argparse
import os
import pandas as pd

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description="Restore SNP-corrected genome based on reference genome. \n \
        Note: SNP-variants show as capitalized charater. \n \
        INDEL-variants show as '?' (insertion) or '_' (deletion). \n \
        e.g. correct_genome.py -R <reference.fa> -vcf_snp <snp.vcf> -vcf_indel <indel.vcf> -O <genome_output> -pre <log_prefix> -sample <sample> \n" 
        )
    parser.add_argument(
        "-R", "--reference", help="reference genome in fasta format", required=True
        )
    parser.add_argument(
        "-gff", "--gff", help="gff file to parse CDS region"
    )
    parser.add_argument(
        "-vcf_snp", "--vcf_snp", help="vcf file with filtered SNP variants",
        )
    parser.add_argument(
        "-snp_restore", "--snp_restore", default = True, help = "if restore snp-corrected genome. default: Yes"
    )
    parser.add_argument(
        "-indel_mark", "--indel_mark", default = False, help = "if mark indel in snp-corrected genome. default: Yes"
    )
    parser.add_argument(
        "-indel_restore", "--indel_restore", default = False, 
        help = "if restore indel in indel-marked genome. default: No"
    )
    parser.add_argument(
        "-vcf_indel", "--vcf_indel", help="vcf file with filtered indel variantes",
    )
    parser.add_argument(
        "-O", "--output", help="output name for the corrected fasta genome", default="restore_genome"
        )
    parser.add_argument(
        "-pre", "--prefix", help="the name of prefix of output file", required=True
    )
    parser.add_argument(
        "-dir", "--directory", help="the directory where output file should be", default = "restore_genome"
    )
    parser.add_argument(
        "-sample", "--sample", required = True, help="the sample name for this restoration"
    )
    args=parser.parse_args()
    return args

def build_dir(dir_name):
    '''Create output directory, write files into'''
    if dir_name not in os.listdir():
        os.mkdir(dir_name)
    return None

def fasta_dict(ref_fasta):
    '''parser reference genome and build dictionary'''
    seq_parse=SeqIO.parse(ref_fasta, "fasta")
    seq_dict=SeqIO.to_dict(seq_parse)
    seq_dict_mut={}
    for entry in seq_dict:
        seq_dict_mut[seq_dict[entry].id] = str(seq_dict[entry].seq)
    print("###### Building reference genome dictionary successfully! ######")
    return seq_dict_mut

def chrom_set(vcf):
    '''chromosome/contig with snp or indel'''
    pvf=pysam.VariantFile(vcf)
    uniq_chr=set()
    for vf_fetch in pvf.fetch():
        uniq_chr.add(vf_fetch.chrom)
    print("###### Build seq id container with variants successfully! ######")
    return sorted(uniq_chr)

def check_char_inseq(seq):
    '''check original sequence has no special character '_' and '?' '''
    if "_" in seq:
        print("+++++ Error '_' in seq !!! ++++++")
    elif "?" in seq:
        print("+++++ Error '?' in seq !!! ++++++")
    else:
        print("***** sequence ready ******")

def snp_restore(args):
    '''restore snp based on reference'''
    output_dir = args.directory
    output_name = args.output
    pre_name = args.prefix
    sample = args.sample
    print("###### reference processing ...... ")
    ref_dict = fasta_dict(args.reference)
    print("###### snp vcf processing ...... ")
    snp_filter = pd.read_table(args.vcf_snp, header=0, dtype={"chrom": str})
    snp_error = open("snp_error.txt", "w")
    for key in ref_dict:
        seq_key = list(ref_dict[key].lower())
        snp_filter_key = snp_filter["chrom" == key]
        for line in snp_filter_key:
            pos = line.split("\t")[1]
            loc = int(pos)-1
            ref = line.split("\t")[2]
            alt = line.split("\t")[3]
            if ref.lower() == seq_key[loc]:
                seq_key[loc] = alt.upper()
            else:
                snp_error.write("vcf and sequence not match -- " + key, pos + "\n")
        sequence = ''.join(seq_key)
        ref_dict[key] = sequence
    
    print("corrected-SNP genome writing ...... ")
    corrected_snp = open(output_dir + "/"+ output_name + "_snp.fa", "w")
    for key in ref_dict:
        corrected_snp.write(">" + key + " dna: " + sample + " restore snp variants from B73_v4 reference genome\n" + ref_dict[key] + "\n")
    corrected_snp.close()
    snp_error.close()
    
def mark_indel(args):
    '''mark indel into genome'''
    output_dir = args.directory
    output_name = args.output
    pre_name = args.prefix
    sample = args.sample
    print("###### read into index file ......")
    indel_filter = pd.read_table(args.vcf_indel, header=0, dtype={"chrom": str})
    indel_mark = open(output_dir + "/indel_mark.log", "w")
    indel_mark.write("chrom\tpos\tref\tmark\talt\t\n")
    indel_error = open("indel_error.txt", "w")
    snp_correct_genome = fasta_dict(output_dir + "/"+ output_name + "_snp.fa")

    for key in snp_correct_genome:
        seq_key = list(snp_correct_genome[key])
        indel_filter_key = indel_filter["chrom" == key]
        for line in indel_filter_key:
            pos = line.split("\t")[1]
            loc = int(pos)-1
            ref = line.split("\t")[2]
            alt = line.split("\t")[3]
            if len(ref) > len(alt):
                '''This is a deletion'''
                ref_len = len(ref)
                stop_loc = loc + ref_len
                fa_seq = ''.join(seq_key[loc:stop_loc])
                ref_seq = ref.lower()
                if fa_seq == ref_seq:
                    X_replace = list("_"*(ref_len-1))
                    seq_key[int(pos):stop_loc] = X_replace
                    repl = seq_key[loc:stop_loc]
                    repl = ''.join(repl)
                    indel_mark.write(chrom + "\t" + pos + "\t" + ref + "\t" + str(repl) + "\t" + alt + "\t\n")
                else:
                    print("++++++ Conflict between reference and sequence! ++++++")
                    indel_error.write(chrom + "\t" + pos + "\t" + ref + "\t" + str(''.join(seq_key[loc:stop_loc])) + "\t\n")

            elif len(ref) < len(alt):
                '''This is an insertion'''
                if seq_key[loc] == ref.lower():
                    '''Check if the ref loc in vcf is the same as the fasta seq'''
                    seq_key[loc] = "?"
                    indel_mark.write(chrom + "\t" + pos + "\t" + ref + "\t?\t" + alt + "\t\n")
                else:
                    print("++++++ Conflict between vcf reference and fasta sequence ! ++++++")
                    indel_error.write(chrom + "\t" + pos + "\t" + ref + "\t" + str(seq_key[loc]) + "\t\n")
            else:
                indel_error.write("This is not a indel")
        sequence = ''.join(seq_key)
        snp_correct_genome[key] = sequence

    print("Marking indel ...... ")
    indel_mark_fa = open(output_dir + "/"+ output_name + "_mark.fa", "w")
    for key in snp_correct_genome:
        indel_mark_fa.write(">" + key + " dna: " + sample + " mark indel variants from B73_v4 reference genome\n" + snp_correct_genome[key] + "\n")
    indel_mark.close()
    indel_error.close()
    indel_mark_fa.close()
    
def parse_gff(args):
    gff_f = args.gff

def indel_restore(args):
    output_dir = args.directory
    output_name = args.output
    pre_name = args.prefix
    sample = args.sample

    snp_indel_mark_fa = output_dir + "/" + output_name + "_indel_mark.fa"
    ref_dict_mark = fasta_dict(snp_indel_mark_fa)
    corrected_indel = open(output_dir + "/" + output_name + "_snp_indel.fa", "w")
    ins_name = output_dir + "/" + pre_name + "_ins.log"
    ins_table = pd.read_table(ins_name, dtype = {"chrom": str}, sep = "\t", header = 0, index_col = 0)
    
    for key in ref_dict_mark:
        corrected_chrom = open(output_dir + "/chr_" + key + "_snp_indel.fa", "w")
        sequence_list = list(ref_dict_mark[key])
        sequence_str = ref_dict_mark[key]
        ins_count = sequence_str.count("?")
        print("? number " + str(ins_count))
        ins_row = len(ins_table[ins_table["chrom"] == key])
        print("row number " + str(ins_row))
        
        i = 0
        for alter_char in ins_table[ins_table["chrom"] == key]["alt"]:
            i += 1
            print(i)
            mark_index = sequence_list.index("?")
            sequence_list[mark_index] = alter_char
        print("###### " + key + " insertion complete ######")
        correct_ins_seq = ''.join(sequence_list)
        correct_indel_seq = correct_ins_seq.replace("_", "")
        print("###### " + key + " deletion complete ######")
        check_char_inseq(correct_indel_seq)
        ref_dict_mark[key] = correct_indel_seq
        corrected_chrom.write(">" + key + " dna: " + sample + " snp-indel-corrected genome from B73_v4 reference genome\n" + ref_dict[key] + "\n")

    print("INDEL-corrected genome writing ...... ")
    for key in ref_dict_mark:
        corrected_indel.write(">" + key + " dna: " + sample + " snp-indel-corrected genome from B73_v4 reference genome\n" + ref_dict[key] + "\n")
    print("Finishing files writing ...... ")
    corrected_indel.close()

def main():
    args = args_parser()
    build_dir(args.directory)
    if args.snp_restore == True:
        snp_restore(args)
    if args.indel_mark == True:
        mark_indel(args)
    if args.indel_restore == True:
        indel_restore(args)
    
#################
#### Run it #####
#################

if __name__ == "__main__":
    main()