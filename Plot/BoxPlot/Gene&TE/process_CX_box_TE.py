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


def Make_index_dic (infile):

    f = open (infile, 'r')
    ls = f.readlines()
    f.close()

    index_dic = {}
    family_dic = {}
    
    for l in ls:
        index_lst = []
        
        es = l.strip().split('\t')
        Chr = es[0]
        family = es[2]
        start = int(es[3])
        stop  = int(es[4])+1
        attr = es[8].split(";")
        ID = attr[0][3:]

        index_lst = [ID, start, stop]
        family_dic [ID] = family

        if Chr not in index_dic:
            index_dic[Chr] = [index_lst]
        else:
            index_dic[Chr].append (index_lst)

    return index_dic, family_dic


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
                ID, start, stop = index

                if (Pos >= start) and (Pos < stop):
                    filtered += 1

                    if ID not in output_dic:
                        output_dic[ID] = [0,0,0,0,0,0]
                    
                    if Type == "CG":
                        output_dic[ID][0] += mC
                        output_dic[ID][1] += C
                    elif Type == "CHG":
                        output_dic[ID][2] += mC
                        output_dic[ID][3] += C
                    elif Type == "CHH":
                        output_dic[ID][4] += mC
                        output_dic[ID][5] += C

        except:
            pass

    f.close()

    return output_dic


def Print_output (output_dic, family_dic, outfile_head):

    o1 = open (f"{outfile_head}_CpG.txt", 'w')
    o2 = open (f"{outfile_head}_CHG.txt", 'w')
    o3 = open (f"{outfile_head}_CHH.txt", 'w')

    for ID in output_dic:

        temp = output_dic[ID]
        family = family_dic[ID]
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

        print(f"{ID}\t{family}\t{CpG_mC}\t{CpG_C}\t{CpG_perc}", file = o1)
        print(f"{ID}\t{family}\t{CHG_mC}\t{CHG_C}\t{CHG_perc}", file = o2)
        print(f"{ID}\t{family}\t{CHH_mC}\t{CHH_C}\t{CHH_perc}", file = o3)

    o1.close()
    o2.close()
    o3.close()


def main():

    f_gff3 = sys.argv[1]            #Bed file - EDTA.TEanno.split.gff3
    f_CX = sys.argv[2]              #Filename - CX report
    outfile_head = sys.argv[3]

    index_dic, family_dic = Make_index_dic (f_gff3)
    output_dic = Process_CX_report (f_CX, index_dic)
    Print_output (output_dic, family_dic, f"{outfile_head}")


if __name__ == "__main__":
    main()
