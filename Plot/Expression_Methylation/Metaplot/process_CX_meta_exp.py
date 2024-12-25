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


def MakeGeneDic (f_gene):

    f = open (f_gene, 'r')
    lst = f.readlines()
    f.close()

    gene_dic = {}
    for l in lst:
        es = l.strip().split('\t')
        fam = es[0]
        gene = es[1]

        gene_dic [gene] = fam
        
    return gene_dic


def MakeFamDic (f_loc, f_gene):

    gene_dic = MakeGeneDic (f_gene)
    
    f = open (f_loc, 'r')
    f.readline() #Skip header
    
    Len = Count_file_line (f_loc) - 1
    chr_lst = [f"Chr{i}" for i in range (1, 13)] #analyzing genes that are located on 12 chromosomes in rice

    fam_dic = {}
    
    for i in range(Len):
        l = f.readline().strip()
        chr, gene, _,_,_,_, is_TE, _, is_representative, *_ = l.split('\t')

        if chr not in chr_lst:
            continue
            
        if gene not in gene_dic:
            continue
        
        if is_representative == "Y" and is_TE == "N":    
            Fam = gene_dic[gene]
            if Fam in fam_dic:
                fam_dic[Fam].append (l)
            else:
                fam_dic[Fam] = [l]
                
    f.close()

    return fam_dic


def printFamDic (Fam_dic, outfile_head):

    o = open (f"{outfile_head}_gene.txt", "w")
    for Fam in Fam_dic:
        Fam_lst = Fam_dic[Fam]
        for l in Fam_lst:
            gene = l.strip().split('\t')[1]
            print (f"{Fam}\t{gene}", file = o)
    o.close()


def Make_index_dic (lst, Updown, Updown_bin, Gene_propbin):

    index_dic = {}
    
    for l in lst:
        index_lst = []
        
        es = l.strip().split('\t')
        Chr = es[0]
        start = int(es[3])
        stop   = int(es[4])
        strand = es[5]
        
        if strand not in ("+","-"):
            continue

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


def Make_FamIndexDic (FamDic, Updown, Updown_bin, Gene_propbin):

    FamIndexDic = {}

    for Fam in FamDic:

        lst = FamDic[Fam]
        IndexDic = Make_index_dic (lst, Updown, Updown_bin, Gene_propbin)
        FamIndexDic[Fam] = IndexDic

    return FamIndexDic


def Initialize_lst (updown, updown_bin, propbin):

    lst = []

    group_no = updown//updown_bin
    total_no = group_no * 2 + propbin

    for i in range(total_no):
        lst.append ([0,0,0,0,0,0])
        
    return lst


def Filter_CX_report (infile, index_dic, output_lst, total_bin):

    f = open (infile, 'r')
    Len = Count_file_line (infile)

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
                if (Pos >= index_lst[0]) and (Pos <= index_lst[-1]):
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


def Print_output (output_lst, outfile_head):

    o1 = open (f"{outfile_head}_CpG.txt", 'w')
    o2 = open (f"{outfile_head}_CHG.txt", 'w')
    o3 = open (f"{outfile_head}_CHH.txt", 'w')

    for i in range(len(output_lst)):

        temp = output_lst[i]
        index = i+0.5
        CpG_mC, CpG_C, CHG_mC, CHG_C, CHH_mC, CHH_C = temp
        if CpG_mC+CpG_C == 0:
            CpG_perc = "NA"
        else:
            CpG_perc = round(Decimal(CpG_mC/(CpG_mC+CpG_C) * 100),2)
    
        if CHG_mC+CHG_C == 0:
            CHG_perc = "NA"
        else:
            CHG_perc = round(Decimal(CHG_mC/(CHG_mC+CHG_C) * 100),2)

        if CHH_mC+CHH_C == 0:
            CHH_perc = "NA"
        else:
            CHH_perc = round(Decimal(CHH_mC/(CHH_mC+CHH_C) * 100),2)

        print(f"{index}\t{CpG_mC}\t{CpG_C}\t{CpG_perc}", file = o1)
        print(f"{index}\t{CHG_mC}\t{CHG_C}\t{CHG_perc}", file = o2)
        print(f"{index}\t{CHH_mC}\t{CHH_C}\t{CHH_perc}", file = o3)

    o1.close()
    o2.close()
    o3.close()


def Process_FamIndexDic (Fam_IndexDic, UpDown, UpDown_bin, Body_propbin, f_CXreport, o_header):

    total_bin = int(UpDown/UpDown_bin)*2 + Body_propbin

    for Fam in Fam_IndexDic:

        IndexDic = Fam_IndexDic[Fam]
        Fam_output_lst = Initialize_lst (UpDown, UpDown_bin, Body_propbin)
        Fam_output_lst = Filter_CX_report (f_CXreport, IndexDic, Fam_output_lst, total_bin)
        Print_output (Fam_output_lst, f"{o_header}_{Fam}")


def main():

    infile_locus = sys.argv[1]       #all locus information (Gene, TE) - all.locus_brief_info.7.0
    infile_gene = sys.argv[2]        #the list of the genes - FG1_exp.txt, ...
    CX_report_file = sys.argv[3]
    outfile_head = sys.argv[4]
    UpDown = int(sys.argv[5])        #2000-bp upstream & downstream regions
    UpDown_bin = int(sys.argv[6])    #100-bp bin for upstream & downstream regions
    Body_propbin = int(sys.argv[7])  #20-proportional bin for Gene body and TE body

    Fam_dic = MakeFamDic (infile_locus, infile_gene)
    printFamDic (Fam_dic, outfile_head)

    Fam_index_dic = Make_FamIndexDic (Fam_dic, UpDown, UpDown_bin, Body_propbin)
    Process_FamIndexDic (Fam_index_dic, UpDown, UpDown_bin, Body_propbin, CX_report_file, outfile_head)

    
if __name__ == "__main__":
    main()

    

    
