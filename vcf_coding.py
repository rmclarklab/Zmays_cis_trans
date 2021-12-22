#! /usr/bin/env python

"""
This script is for parsing variants on gene coding region by providing filtered 
vcf file and gene feature file (gff).
NOTE: the vcf file must have been filtered.
By Meiyuan Ji, meiyuan.ji@utah.edu
Last update 11/19/2019.
"""

import pysam
import argparse
import subprocess
import os

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="filtered SNPs on exons with other gene information. \n\
        Usage: vcf_coding.py -vcf <vcf_file> -gtf <sorted_gtf> -O <ouput_name>")
    parser.add_argument("-vcf", "--vcf", help="vcf file with passed snp. ")
    parser.add_argument("-gtf", "--gtf", help="sorted gtf file with tabix index. ")
    parser.add_argument("-O", "--output", help="output vcf file with variants on exon region. ")
    args = parser.parse_args()
    return args

def gtf_contig(args):
    '''gtf chromosome/contig '''
    gtf = args.gtf
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    contig_set = set()
    gtf_ps = gtf_handle.fetch()
    for ginfo in gtf_ps:
        contig = ginfo.contig
        contig_set.add(contig)
    return sorted(contig_set)

def chrom_set(args):
    '''vcf chromosome/contig with snp'''
    vcf = args.vcf
    pvf=pysam.VariantFile(vcf)
    uniq_chr=set()
    for vf_fetch in pvf.fetch():
        uniq_chr.add(vf_fetch.chrom)
    print("###### Build seq id container with variants successfully! ######")
    return sorted(uniq_chr)

def gene_pos(args):
    '''acquire all protein coding genes position information in the gtf file. '''
    print("###### Processing gene position information ...... ")
    # gene_dict_file = open("gene_dict.txt", "w")
    gene_t = ['pseudogene','protein_coding']
    gene_dict = {}
    contig_set = gtf_contig(args)
    gtf = args.gtf
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    for contig in contig_set:
        gene_dict[contig] = {}
        for gr in gtf_handle.fetch(contig):
            if (gr.feature == "gene" and gr.gene_biotype in gene_t):
                # positions in pysam are all 0-based
                # start is the gtf_start-1, and end is the value gtf_end (since the nature of python as right excluded)
                gene_dict[contig][gr.gene_id] = [gr.start, gr.end]
    # gene_dict_file.write(str(gene_dict))
    return gene_dict

def exon_dict(args, gene_dict):
    '''exon dictionary with all transcript id each gene_id information in the gtf file. '''
    print("###### Processing transcript id of each gene ......")
    # exon_dict = open("exon_dict.txt", "w")
    gtf = args.gtf
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    exon_dict = {}
    for contig in gene_dict:
        exon_dict[contig] = {}
        for gid in gene_dict[contig]:
            # gene start position is 0-based
            exon_dict[contig][gid] = []
            start = gene_dict[contig][gid][0]
            end = gene_dict[contig][gid][-1]
            for gexon in gtf_handle.fetch(contig, start, end):
                # find all exons 
                if gexon.feature == "exon":
                    # exon_set.add(gexon.transcript_id)
                    exon_dict[contig][gid].append([gexon.start, gexon.end])
    # the exon_dict includes keys exon_dict[contig][gene_id] = [[start1, end1], [start2, end2], ......]
    # exon_dict.write(str(exon_dict))
    return exon_dict

def exon_pos(args, gene_dict, exon_dict):
    '''exon position for each gene id in the gtf file. '''
    print("###### Processing exon position of each gene ...... ")
    gtf = args.gtf
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    exon_pos_dict = {}
    for contig in exon_dict:
        exon_pos_dict[contig] = {}
        for gid in exon_dict[contig]:
            # here get the gene positions #
            start = gene_dict[contig][gid][0]
            end = gene_dict[contig][gid][-1]
            # here get the gene positions # 
            exon_pos_dict[contig][gid] = {}
            for tid in exon_dict[contig][gid]:
                exon_pos_dict[contig][gid][tid] = []
                for exonid in gtf_handle.fetch(contig, start, end):
                    if (exonid.feature == 'exon' and exonid.transcript_id == tid):
                        # here, the start are all 0-based 
                        # end are 1-based because the the exclusive nature of python
                        exon_pos_dict[contig][gid][tid].append([exonid.start+1, exonid.end])
    return exon_pos_dict

def parse_pos_on_exon(pos, exon_dict, contig, gid):
    # print("###### Processing variants position ...... ")
    for pos_l in exon_dict[contig][gid]:
        # isoforms_list = exon_pos_dict[contig][gid][tid]
        # for iso in isoforms_list:
        if pos_l[0] <= pos <= pos_l[-1]:
            return True
            break

def coding_v(args, gene_dict, exon_dict):
    # parser gtf file 
    gtf = args.gtf
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    gtf_chr = gtf_contig(args)
    # sorted chromosome set with passed variants on it 
    vcf_contig_set = chrom_set(args)
    # parser vcf file
    vcf = args.vcf
    ###
    vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    vcf_w = open(output + ".vcf", "w")
    vcf_w.write(hd)
    ### 

    for contig in vcf_contig_set:
        # contig means contig
        print("processing contig " + contig)
        if contig in gtf_chr:
            for gid in gene_dict[contig]:
                var_dict = {}
                vpos = []
                vref = []
                valt = []
                # the positions are 0-based, used for gtf fetch
                start = gene_dict[contig][gid][0]
                end = gene_dict[contig][gid][-1]
                # get back to the original value as in gtf file
                gstart = start+1
                gend = end
                gnpos = str(gstart) + "-" + str(gend)
                # fetch are 0-based
                # and the positions are originally 0-based, no need to change 
                for var in vcf_handle.fetch(contig, start, end):
                    # the position in vcf are 1-based, but I read it with fetch method is 0-based
                    if parse_pos_on_exon(var.pos + 1, exon_dict, contig, gid):
                        var_dict[var.pos+1] = [var.ref, var.alt]
                        vpos.append(var.pos+1)
                        vref.append(var.ref)
                        valt.append(var.alt)
                        
                        ###
                        vcf_w.write(str(var) + "\n")
    vcf_w.close()

def sort_vcf(args):
    outputl = args.output
    output = outputl + ".vcf"
    sort_command = f"bcftools sort -o {output}.gz -O z {output}"
    subprocess.call(sort_command, shell = True)
    index_command = f"bcftools index -t --threads 4 {output}.gz"
    subprocess.call(index_command, shell = True)
    gzvcf_file = output + ".gz"
    if gzvcf_file in os.listdir():
        rm_command = f"rm {output}"
        subprocess.call(rm_command, shell = True)
        print("Deleted the uncompressed vcf file. ")

def main():
    args = args_parser()
    # get gene dictionary as gene_dict[contig][gene_id] = [start, end]
    gene_dict = gene_pos(args)
    # get exon dictionary as exon_dict[contig][gene_id] = [transcript_ids]
    exon_d = exon_dict(args, gene_dict)
    # get exon position dictionary as exon_pos[contig][gid][tid] = [[start1, end1], [start2, end2], .......]
    # exon_p = exon_pos(args, gene_dict, exon_d)
    # write the variants on exon regions as dictionary to output file
    coding_v(args, gene_dict, exon_d)
    sort_vcf(args)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()