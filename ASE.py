"""Update 02/18/2020
Modified for filtering genes with bad SNP-position, 
and count reads based on individual SNP sites. """

"""This script is to filter bads genes for each genotype 
according to the read count results. I filter out bad genes 
based on 2 conditions: 
1). mis-mapped reads which greater than 5 will be list as candidate bad genes;
2). mis-mapped reads account for some proportion (cutoff) will at last list as bad genes. """

import pandas as pd
import argparse
import os
import pysam 

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Perform allele-specific read count on gene-basis. \n\
        Usage: python ASE.py ")
    parser.add_argument("-v", "--vcf", help = "VCF file with only variants information on exon regions. ")
    parser.add_argument("-gtf", "--gtf", help = "gene annotation file in gtf format. ")
    parser.add_argument("-bam", "--bam", help = "the bam file for read counts. ")
    parser.add_argument("-O", "--output", help = "output name ")
    args = parser.parse_args()
    return args

def gene_info(args):
    '''For each gene id, get the gene oritation from gtf file. Return as a dictionary. '''
    ginfo = {}
    gtf_handle = pysam.TabixFile(args.gtf, parser = pysam.asGTF())
    ### only consider genes with variants
    for gene in gtf_handle.fetch():
        if (gene.feature == "gene") and (gene.gene_biotype in ["pseudogene", "protein_coding"]):
            gtig = gene.contig
            gid = gene.gene_id
            gstrand = gene.strand
            gstart = gene.start
            gend = gene.end
            ginfo[gid] = [gstrand, gtig, gstart, gend]
    return ginfo

def read_strand(read, gene_id, strand):
    if strand == "+":
        if ((read.is_read1 and read.is_reverse) or (read.is_read2 and read.is_reverse == False)):
            return True
        else:
            return False
    elif strand == "-":
        if ((read.is_read1 and read.is_reverse == False) or (read.is_read2 and read.is_reverse)):
            return True
        else:
            return False
    else:
        return False

def allele_gene_count(args, gid_inf):
    '''Read counts based on SNP sites. '''
    ### open vcf, bam and gtf
    vcf_handle = pysam.TabixFile(args.vcf, parser = pysam.asVCF())
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    gtf_handle = pysam.TabixFile(args.gtf, parser = pysam.asGTF())
    ### output file
    filew = open(args.output + ".txt", "w")
    filew.write("gene\ttotal_count\tref_count\talt_count\tnon_count\n")
    ### get gene strand dictionary given gtf file
    ginfo = gene_info(args)
    ### 
    for gg in ginfo: # gid:[gstrand,gtig, gstart,gend]
        gstrand = ginfo[gg][0]
        gtig = ginfo[gg][1]
        gstart = ginfo[gg][2]
        gend = ginfo[gg][3]
        total_read = set()
        ref_read = set()
        alt_read = set()
        non_read = set()
        for var in vcf_handle.fetch(gtig, gstart, gend):
            vpos = var.pos
            ref = var.ref
            alt = var.alt
            ### based on SNP positions
            for read_column in bam_handle.pileup(gtig, vpos-1, vpos):
                if read_column.reference_pos == vpos-1: # important, cannot delete
                    for bam_read in read_column.pileups:
                        read_name = bam_read.alignment.query_name
                        # uniquely mapped, mapped to exon, read strand
                        if (("NH", 1) in bam_read.alignment.tags and bam_read.query_position != None and read_strand(bam_read.alignment, gg, gstrand)):
                            ###
                            total_read.add(read_name)
                            ###
                            allele = bam_read.alignment.query_sequence[bam_read.query_position]
                            # print("chromosome " + chr  + " position " + str(p) + " ref " + ref + " alt " + alt + " allele " + allele)
                            if allele == ref:
                                ref_read.add(read_name)
                            elif allele == alt:
                                alt_read.add(read_name)
                            else:
                                non_read.add(read_name)
            ###
            total_read_count = len(total_read)
            ref_read_count = len(ref_read)
            alt_read_count = len(alt_read)
            non_read_count = len(non_read)
            filew.write(f"{gg}\t{total_read_count}\t{ref_read_count}\t {alt_read_count}\t{non_read_count}\n")
    filew.close()

def main():
    args = args_parser()
    gid_inf = gene_info(args)
    allele_gene_count(args, gid_inf)

################
#### Run it ####
################

if __name__ == "__main__":
    main()