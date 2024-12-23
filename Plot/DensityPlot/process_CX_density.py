import sys
from decimal import *

def Count_file_line (infile):

    i = 0
    f = open (infile, 'r')
    for i, line in enumerate(f):
      pass
    Len = i+1
    
    f.close()

    return Len


def CX_to_dic (f_CX, dic):

    f = open (f_CX, 'r')
    Len = Count_file_line (f_CX)

    for i in range (Len):
        l = f.readline()
        es = l.strip().split('\t')

        try:
            chr = int(es[0][3:])
            pos  = int(es[1])
            mC   = int(es[3])
            C    = int(es[4])
            Context = es[5]

            if chr not in dic:
                dic[chr] = {}
            if pos not in dic[chr]:
                dic[chr][pos] = [0,0,0,0,0,0]

            if Context == "CG":
                dic[chr][pos][0] += mC
                dic[chr][pos][1] += C
            elif Context == "CHG":
                dic[chr][pos][2] += mC
                dic[chr][pos][3] += C
            elif Context == "CHH":
                dic[chr][pos][4] += mC
                dic[chr][pos][5] += C
        
        except:
            pass

    f.close()

    return dic


def Print_bed_w1 (dic, f_o_header):

    o_CG  = open (f"{f_o_header}_CG_w1.bed", 'w')
    o_CHG = open (f"{f_o_header}_CHG_w1.bed", 'w')
    o_CHH = open (f"{f_o_header}_CHH_w1.bed", 'w')

    for chr in dic:
        pos_dic = dic[chr]

        for pos in pos_dic:
            start = pos - 1
            end   = pos
            CG_mC, CG_C, CHG_mC, CHG_C, CHH_mC, CHH_C = pos_dic[pos]

            if CG_mC + CG_C > 0:
                CG_perc = round(Decimal(CG_mC/(CG_mC+CG_C) * 100),2)
                print (f"{chr}\t{start}\t{end}\t{CG_mC}\t{CG_C}\t{CG_perc}", file = o_CG)

            if CHG_mC + CHG_C > 0:
                CHG_perc = round(Decimal(CHG_mC/(CHG_mC+CHG_C) * 100),2)
                print (f"{chr}\t{start}\t{end}\t{CHG_mC}\t{CHG_C}\t{CHG_perc}", file = o_CHG)

            if CHH_mC + CHH_C > 0:
                CHH_perc = round(Decimal(CHH_mC/(CHH_mC+CHH_C) * 100),2)
                print (f"{chr}\t{start}\t{end}\t{CHH_mC}\t{CHH_C}\t{CHH_perc}", file = o_CHH)

    o_CG.close()
    o_CHG.close()
    o_CHH.close()
    

def Initialize_chr (f_chr, window):
    
    f = open (f_chr, 'r')
    ls = f.readlines()
    f.close()

    dic_chr = {}
    dic_output = {}
    for l in ls:
        es = l.strip().split('\t')
        chr = int (es[0])
        chr_len = int (es[1])
        dic_chr[chr] = chr_len

        chr_bin = chr_len//window
        temp = [[0,0] for _ in range(chr_bin+1)]
        dic_output[chr] = temp

    return dic_chr, dic_output


def bed_w1_to_w50 (f_o_header, f_chr):

    for Context in ["CG","CHG","CHH"]:
        f_i = f"{f_o_header}_{Context}_w1.bed"
        f_o = f"{f_o_header}_{Context}_w50.bed"
        
        f = open (f_i, 'r')
        o = open (f_o, 'w')
        Len = Count_file_line (f_i)

        dic_chr, dic_output = Initialize_chr (f_chr, 50)
        
        for x in range(Len):
            l = f.readline()
            es = l.strip().split('\t')
            chr   = int(es[0])
            pos   = int(es[1])
            mC    = int(es[3])
            C     = int(es[4])

            ID = pos//50
            dic_output[chr][ID][0] += mC
            dic_output[chr][ID][1] += C
        f.close()

        for chr in dic_output:
            temp = dic_output[chr]
            chr_len = dic_chr[chr]
            chr_bin = len(temp)

            for i in range(chr_bin):
                start = i*50
                end   = (i+1)*50
                if end > chr_len:
                    end = chr_len
                mC    = temp[i][0]
                C     = temp[i][1]
                
                if mC + C > 0:
                    perc = round(Decimal(mC/(mC+C) * 100),2)
                    print (f"{chr}\t{start}\t{end}\t{mC}\t{C}\t{perc}", file = o)
        o.close()


def main():

    #Note that CX report is one-based, while .bed is zero-based

    f_chr = sys.argv[1]              #infile - chromosome
    f_o_header  = sys.argv[2]        #outfile_header
    lst_f_CX = sys.argv[3:]          #Files (CX report) to be merged into one bed file

    dic = {}
    for f_CX in lst_f_CX:
        dic = CX_to_dic (f_CX, dic) 

    Print_bed_w1 (dic, f_o_header)
    bed_w1_to_w50 (f_o_header, f_chr) 


if __name__ == "__main__":
    main()
