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


def Make_Gene_lst (infile):

    f = open (infile, 'r')
    f.readline() #Skip header
    
    Len = Count_file_line (infile) - 1

    Gene_lst = []
    
    for i in range(Len):
        l = f.readline()
        es = l.strip().split('\t')

        is_TE  = es[6]
        is_representative = es[8]

        if (is_representative == "Y") and (is_TE == "N"):
            Gene_lst.append(l)

    return Gene_lst


def Make_index_dic (lst):

    index_dic = {}
    
    for l in lst:
        index_lst = []
        
        es = l.strip().split('\t')
        Chr = es[0]
        locus = es[1]
        start = int(es[3])
        stop  = int(es[4])+1

        if Chr not in ["ChrUn", "ChrSy"]:

            index_lst = [locus, start, stop]

            if Chr not in index_dic:
                index_dic[Chr] = [index_lst]
            else:
                index_dic[Chr].append (index_lst)

    return index_dic


def Process_CX_report (infile, index_dic):

    f = open (infile, 'r')
    Len = Count_file_line (infile)
    
    filtered = 0
    output_dic = {}

    for i in range (Len):
        l = f.readline()
        es = l.strip().split('\t')

        Chr  = es[0]
        Pos  = int(es[1])
        mC   = int(es[3])
        C    = int(es[4])
        Type = es[5]

        try:
            index_lst = index_dic[Chr]
            for index in index_lst:
                locus, start, stop = index

                if (Pos >= start) and (Pos < stop):
                    filtered += 1

                    if locus not in output_dic:
                        output_dic[locus] = [0,0,0,0,0,0]
                    
                    if Type == "CG":
                        output_dic[locus][0] += mC
                        output_dic[locus][1] += C
                    elif Type == "CHG":
                        output_dic[locus][2] += mC
                        output_dic[locus][3] += C
                    elif Type == "CHH":
                        output_dic[locus][4] += mC
                        output_dic[locus][5] += C

        except:
            pass
    
    f.close()

    return output_dic


def Print_output (output_dic, outfile_head):

    o1 = open (f"{outfile_head}_CpG.txt", 'w')
    o2 = open (f"{outfile_head}_CHG.txt", 'w')
    o3 = open (f"{outfile_head}_CHH.txt", 'w')

    for Locus in output_dic:

        temp = output_dic[Locus]
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

        print(f"{Locus}\t{CpG_mC}\t{CpG_C}\t{CpG_perc}", file = o1)
        print(f"{Locus}\t{CHG_mC}\t{CHG_C}\t{CHG_perc}", file = o2)
        print(f"{Locus}\t{CHH_mC}\t{CHH_C}\t{CHH_perc}", file = o3)

    o1.close()
    o2.close()
    o3.close()


def main():

    f_locus = sys.argv[1]           #all locus information (Gene, TE) - all.locus_brief_info.7.0
    f_CX = sys.argv[2]              #Filename - CX report
    outfile_head = sys.argv[3]

    Gene_lst = Make_Gene_lst (f_locus)
    Gene_index_dic = Make_index_dic (Gene_lst)
    Gene_output_dic = Process_CX_report (f_CX, Gene_index_dic)
    Print_output (Gene_output_dic, f"{outfile_head}_Gene")


if __name__ == "__main__":
    main()

    

    
