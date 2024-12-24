import argparse
from decimal import *

def setup_argparse():
    parser = argparse.ArgumentParser()

    parser.add_argument ('region_file', help="bed file to define genomic regions of interest")
    parser.add_argument ('CX_file', help= "CX_report")
    parser.add_argument ('out_file')
    parser.add_argument ('--chr_header', type=str, default="", 
                         help = "header of chromosome that needs to be removed from due to the difference with region_file")
    args = parser.parse_args()

    return args


def make_index_dict (region_file):

    index_dict = {}
    
    with open (region_file, 'r') as f:
        for line in f:
            chr, start, end = line.strip().split('\t')
            start = int(start)
            end = int(end)

            if chr not in index_dict: 
                index_dict[chr] = []
            index_dict[chr].append ([start, end])

    return index_dict


def filter_CX (region_file, CX_file, temp_file, chr_header): 

    indexDict = make_index_dict (region_file)
    o = open (temp_file, 'w')

    with open (CX_file) as f:
        for line in f:
            
            chr, pos, _, mC, C, context, _ = line.strip().split('\t')
            chr = chr.replace(chr_header, "")
            pos = int(pos) - 1 #CX report uses 1-based coordinate
            mC  = int(mC)
            C   = int(C)

            for regionStart, regionEnd in indexDict[chr]:
                if (pos >= regionStart) and (pos < regionEnd):
                    if mC + C > 0:
                        print (f"{chr}\t{regionStart}\t{regionEnd}\t{mC}\t{C}\t{context}", file = o, flush=True)

    o.close()


def make_region_dict (temp_file):

    region_dict = {}

    with open (temp_file) as f:
        for line in f:
            
            chr, start, end, mC, C, context = line.strip().split('\t')
            mC = int(mC)
            C = int(C)

            if chr not in region_dict:
                region_dict[chr] = {}
            if start not in region_dict[chr]:
                region_dict[chr][start] = [end, 
                                          {"CG":  [0,0],
                                           "CHG": [0,0],
                                           "CHH": [0,0]}]

            region_dict[chr][start][1][context][0] += mC
            region_dict[chr][start][1][context][1] += C

    return region_dict


def process_CX (temp_file, out_file):

    region_dict = make_region_dict (temp_file)
    o = open (out_file, 'w')

    for chr in region_dict:
        for start in region_dict[chr]:
            end = region_dict[chr][start][0]
            methyl_dict = region_dict[chr][start][1]

            for context in methyl_dict:
                mC, C = methyl_dict[context]
                if mC+C != 0:
                    mC_rate = round(Decimal(mC/(mC+C)*100),2)
                    print (f"{chr}\t{start}\t{end}\t{mC}\t{C}\t{mC_rate}\t{context}", file = o)

    o.close()


def main():

    args = setup_argparse()
    region_file = args.region_file
    CX_file = args.CX_file
    out_file = args.out_file
    chr_header = args.chr_header

    temp_file = f"{out_file}.temp"

    filter_CX (region_file, CX_file, temp_file, chr_header)
    process_CX (temp_file, out_file)


if __name__ == "__main__":
    main()
