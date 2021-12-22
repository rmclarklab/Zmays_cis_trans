import pandas as pd
import argparse
import os
import pysam 

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="filter bad genes for each genotype. \n\
        Usage: SNP_filter.py -dir <read_count> -is_ref -O ")
    parser.add_argument("-dir", "--directory", help = "the directory with allelic read count of given genotype. ")
    parser.add_argument("-is_ref", "--is_ref", help = "assign the genotype as reference genome. ", action = "store_true", default = False)
    parser.add_argument("-O", "--output", help = "output name ")
    parser.add_argument("-cutoff", "--cutoff", default = 5, help = "the percentage cutoff value to filter snp sites. ", type = int)
    args = parser.parse_args()
    return args

def count_sum(args):
    dir = args.directory
    fs = [f for f in os.listdir(dir) if f.endswith(".txt")]
    f_total = []
    f_ref = []
    f_alt = []
    f_non = []
    for f in fs:
        snp_allele = pd.read_table(os.path.join(dir, f), sep = "\t", header = 0, index_col = ["chrom", "variants_pos"])
        snp_total = snp_allele[["total_count"]]
        snp_ref = snp_allele[["ref_count"]]
        snp_alt = snp_allele[["alt_count"]]
        snp_non = snp_allele[["non_count"]]
        f_total.append(snp_total)
        f_ref.append(snp_ref)
        f_alt.append(snp_alt)
        f_non.append(snp_non)
    total_df = pd.concat(f_total, axis = 1)
    ref_df = pd.concat(f_ref, axis = 1)
    alt_df = pd.concat(f_alt, axis = 1)
    non_df = pd.concat(f_non, axis = 1)
    ###
    total_sum_col = total_df.sum(axis = 1)
    ref_sum_col = ref_df.sum(axis = 1)
    alt_sum_col = alt_df.sum(axis = 1)
    non_sum_col = non_df.sum(axis = 1)
    ####
    df = pd.DataFrame(columns = ["total", "ref", "alt", "non"])
    df["total"] = total_sum_col
    df["ref"] = ref_sum_col
    df["alt"] = alt_sum_col
    df["non"] = non_sum_col
    return df

def bad_snp(args, df):
    '''Get mis-ratio of filtered genes against total-count. '''
    is_ref = args.is_ref
    cutoff = args.cutoff
    output = args.output
    #
    if is_ref == True:
        df["mis_ratio"] = (df["alt"] + df["non"])/df["total"] * 100
    else:
        df["mis_ratio"] = (df["ref"] + df["non"])/df["total"] * 100
    ########### cutoff for potential problematic genes set here ##############
    bad_snp = df[df.mis_ratio > cutoff]
    bad_snp.to_csv(output+".txt", sep = "\t", index = True)
    return bad_gid

def main():
    args = args_parser()
    df_sum = count_sum(args)
    bad_snp(args, df_sum)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
