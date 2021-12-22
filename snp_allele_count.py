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
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="filter bad genes for each genotype. \n\
        Usage: SNP_count.py -dir <read_count> -gt <genotype> -is_ref")
    parser.add_argument("-vcf", "--vcf", help = "VCF file with only variants information on exon regions. ")
    parser.add_argument("-bam", "--bam", help = "the bam file for read counts. ")
    parser.add_argument("-gtf", "--gtf", help = "gene annotation file in gtf format. ")
    parser.add_argument("-O", "--output", help = "output name ")
    args = parser.parse_args()
    return args

def gene_strand(args):
    '''For each gene id, get the gene oritation from gtf file. Return as a dictionary. '''
    gstrand_dict = {}
    gtf_handle = pysam.TabixFile(args.gtf, parser = pysam.asGTF())
    ### only consider genes with variants
    for gene in gtf_handle.fetch():
        if (gene.feature == "gene") and (gene.gene_biotype in ["pseudogene", "protein_coding"]):
            gid = gene.gene_id
            gstrand = gene.strand
            gstrand_dict[gid] = gstrand
    return gstrand_dict

def read_strand(read, gstrand, gene_id):
    if gstrand[gene_id] == "+":
        if ((read.is_read1 and read.is_reverse) or (read.is_read2 and read.is_reverse == False)):
            return True
        else:
            return False
    elif gstrand[gene_id] == "-":
        if ((read.is_read1 and read.is_reverse == False) or (read.is_read2 and read.is_reverse)):
            return True
        else:
            return False
    else:
        return False

def snp_count(args):
    '''Read counts based on SNP sites. '''
    ### open vcf, bam and gtf
    vcf_handle = pysam.TabixFile(args.vcf, parser = pysam.asVCF())
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    ### output file
    snp_file = open(args.output + ".txt", "w")
    snp_file.write("chrom\tvariants_pos\ttotal_count\tref_count\talt_count\tnon_count\n")
    ### get gene strand dictionary given gtf file
    ### 
    gstrand = gene_strand(args)
    for var in vcf_handle.fetch():
        vcontig = var.contig
        vpos = var.pos # 0-based
        ref = var.ref
        alt = var.alt
        ### based on SNP positions
        for read_column in bam_handle.pileup(vcontig, vpos-1, vpos):
            if read_column.reference_pos == p-1: # important, cannot delete
                for bam_read in read_column.pileups:
                    read_name = bam_read.alignment.query_name
                    # uniquely mapped, mapped to exon, read strand
                    if (("NH", 1) in bam_read.alignment.tags and bam_read.query_position != None and read_strand(bam_read.alignment, gstrand, gid)):
                        ###
                        total_read.add(read_name)
                        ###
                        allele = bam_read.alignment.query_sequence[bam_read.query_position]
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
            snp_file.write("{0:s}\t{1:s}\t{2:d}\t{3:d}\t{4:d}\t{5:d}\t{6:d}\n".format(chr, gid, p, total_read_count, ref_read_count, alt_read_count, non_read_count))
    snp_file.close()

def main():
    args = args_parser()
    snp_count(args)

################
#### Run it ####
################

if __name__ == "__main__":
    main()