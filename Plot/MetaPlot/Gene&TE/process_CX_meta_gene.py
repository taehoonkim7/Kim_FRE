import sys
from decimal import *

def count_file_line (infile):

    i = 0
    f = open (infile, 'r')
    for i, line in enumerate(f):
      pass
    Len = i+1
    
    f.close()

    return Len


def make_gene_lst (infile):

    f = open (infile, 'r')
    f.readline() #Skip header
    
    Len = count_file_line (infile) - 1

    gene_lst = []
    
    for i in range(Len):
        l = f.readline()
        es = l.strip().split('\t')

        is_TE  = es[6]
        is_representative = es[8]

        if (is_representative == "Y") and (is_TE == "N"):
            gene_lst.append(l)

    return gene_lst
 

def initialize_lst (updown, updown_bin, propbin):

    lst = []

    group_no = updown//updown_bin
    total_no = group_no * 2 + propbin

    for i in range(total_no):
        lst.append ([0,0,0,0,0,0])
        
    return lst


def make_index_dic (lst, Updown, Updown_bin, Gene_propbin):

    index_dic = {}
    
    for l in lst:
        index_lst = []
        
        es = l.strip().split('\t')
        Chr = es[0]
        start = int(es[3])
        stop  = int(es[4])+1
        strand = es[5]

        gap = stop - start
        
        Up_start = start - Updown 
        Down_end = stop + Updown
        bin = Updown//Updown_bin

        for i in range(bin):
            index = Up_start + Updown_bin * i
            index_lst.append (index)
        for i in range(Gene_propbin):
            index = start + int(gap * (i/Gene_propbin))
            index_lst.append (index)
        for i in range(bin+1):
            index = stop + Updown_bin * i
            index_lst.append (index)

        if Chr not in index_dic:
            index_dic[Chr] = [[strand,index_lst]]
        else:
            index_dic[Chr].append ([strand,index_lst])

    return index_dic


def filter_CX_report (infile, index_dic, output_lst, total_bin):

    f = open (infile, 'r')
    Len = count_file_line (infile)

    filtered = 0
    for i in range (Len):
        l = f.readline()
        es = l.strip().split('\t')

        Chr  = es[0]
        Pos  = int(es[1])
        mC   = int(es[3])
        C    = int(es[4])
        Type = es[5]

        try:
            index_lst_lst = index_dic[Chr]
            for strand, index_lst in index_lst_lst:
                if (Pos >= index_lst[0]) and (Pos < index_lst[-1]):
                    filtered += 1
                    
                    for j in range(len(output_lst)):    
                        if (Pos >= index_lst[j]) and (Pos < index_lst[j+1]):
                            if strand == "+":
                                index = j
                            elif strand == "-":
                                index = total_bin-1-j
                            
                            if Type == "CG":
                                output_lst[index][0] += mC
                                output_lst[index][1] += C
                            elif Type == "CHG":
                                output_lst[index][2] += mC
                                output_lst[index][3] += C
                            elif Type == "CHH":
                                output_lst[index][4] += mC
                                output_lst[index][5] += C
        except:
            pass

    f.close()

    return output_lst


def print_output (output_lst, outfile_head):

    o1 = open (f"{outfile_head}_CpG.txt", 'w')
    o2 = open (f"{outfile_head}_CHG.txt", 'w')
    o3 = open (f"{outfile_head}_CHH.txt", 'w')

    for i in range(len(output_lst)):

        temp = output_lst[i]
        index = i+0.5
        CpG_mC, CpG_C, CHG_mC, CHG_C, CHH_mC, CHH_C = temp

        if CpG_mC+CpG_C == 0:
            CpG_perc = 0
        else:
            CpG_perc = round(Decimal(CpG_mC/(CpG_mC+CpG_C) * 100),2)

        if CHG_mC+CHG_C == 0:
            CHG_perc = 0
        else:
            CHG_perc = round(Decimal(CHG_mC/(CHG_mC+CHG_C) * 100),2)

        if CHH_mC+CHH_C == 0:
            CHH_perc = 0
        else:
            CHH_perc = round(Decimal(CHH_mC/(CHH_mC+CHH_C) * 100),2)

        print(f"{index}\t{CpG_mC}\t{CpG_C}\t{CpG_perc}", file = o1)
        print(f"{index}\t{CHG_mC}\t{CHG_C}\t{CHG_perc}", file = o2)
        print(f"{index}\t{CHH_mC}\t{CHH_C}\t{CHH_perc}", file = o3)

    o1.close()
    o2.close()
    o3.close()


def main():

    infile = sys.argv[1]             #all locus information (Gene, TE) - all.locus_brief_info.7.0
    CX_report_file = sys.argv[2]
    outfile_head = sys.argv[3]
    UpDown = int(sys.argv[4])        #2000-bp upstream & downstream regions
    UpDown_bin = int(sys.argv[5])    #100-bp bin for upstream & downstream regions
    Body_propbin = int(sys.argv[6])  #20-proportional bin for gene body

    total_bin = int(UpDown/UpDown_bin)*2 + Body_propbin

    gene_lst = make_gene_lst (infile)
    gene_index_dic = make_index_dic (gene_lst, UpDown, UpDown_bin, Body_propbin)
    gene_output_lst = initialize_lst (UpDown, UpDown_bin, Body_propbin)
    gene_output_lst = filter_CX_report (CX_report_file, gene_index_dic, gene_output_lst, total_bin)
    print_output (gene_output_lst, f"{outfile_head}_gene")

main()

    

    

