import pandas as pd
import random
import argparse

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Subsample read count using random generator \n\
        Usage: random_count.py -I <input> -O <output>")
    parser.add_argument("-I", "--input", help = "the input file with read count")
    parser.add_argument("-O", "--output", help = "the output name")
    args = parser.parse_args()
    return args

def random_generator(args):
    infile = args.input
    outfile = args.output
    df = pd.read_table(infile, sep = "\t", header = 0, index_col = 0)
    for row in df.index:
        for col in df.columns:
            df.at[row, col] = rand(df.at[row, col])
    print("###### writing to a new file ...... ")
    df.to_csv(outfile + ".txt", sep = "\t", index = True)

def rand(num):
    new_num = 0
    for i in range(num):
        if random.randrange(2) == 1:
            new_num += 1
    return new_num

def main():
    args = args_parser()
    random_generator(args)

##################
##### Run it #####
##################

if __name__ == "__main__":
    main()