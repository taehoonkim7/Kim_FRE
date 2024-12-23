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


def Filter_CX_report (infile, outfile, cutoff):

    f = open (infile, 'r')
    o = open (outfile, 'w')

    Len = Count_file_line (infile)
    
    for i in range(Len):
        l = f.readline()
        es = l.strip().split('\t')

        mC = int(es[3])
        C = int(es[4])
        total = mC+C

        if total >= cutoff:
            print (l, end = "", file = o)

    f.close()
    o.close()   


def Make_Chr_ID_dic (Chr_file):

    f = open (Chr_file, 'r')
    ls = f.readlines()
    f.close()

    Chr_ID_dic = {}
    for l in ls:
        es = l.strip().split(" ")
        Chr = es[0]
        Len = int(es[1])

        Chr_ID_dic[Chr] = Len

    return Chr_ID_dic
    

def Initialize_Chr_dic (Chr_file, Step):

    Chr_ID_dic = Make_Chr_ID_dic (Chr_file)
    Chr_dic = {}

    for Chr in Chr_ID_dic.keys():
        Len = Chr_ID_dic[Chr]
        loop = Len//Step

        temp_lst = []
        for i in range(loop+1):
            temp_lst.append ([0,0,0,0,0,0])

        Chr_dic[Chr] = temp_lst

    return Chr_dic


def Compress_CX_report (infile, chr_file, window, step):

    Chr_dic = Initialize_Chr_dic (chr_file, step)
    
    f = open (infile, 'r')
    Len = Count_file_line (infile)

    for i in range(Len):
        l = f.readline()
        es = l.strip().split('\t')

        Chr  = es[0]
        Pos  = int(es[1])
        mC   = int(es[3])
        C    = int(es[4])
        Type = es[5]

        end = Pos//step + 1
        gap = window//step
        start = end - gap
        if start < 0:
            start = 0        

        if Chr not in Chr_dic.keys():
            pass
        else:
            Chr_lst = Chr_dic[Chr]

            for i in range(start, end):
                temp_lst = Chr_lst[i]
        
                if Type == "CG":
                    temp_lst [0] += mC
                    temp_lst [1] += C
                elif Type == "CHG":
                    temp_lst [2] += mC
                    temp_lst [3] += C
                elif Type == "CHH":
                    temp_lst [4] += mC
                    temp_lst [5] += C
                else:
                    print (f"Compress_CX_report ERROR: Type = {Type}", flush = True)

    return Chr_dic
        
        
def Print_Chr_dic (outfile_head, Chr_dic, window, step):

    o1 = open (f"{outfile_head}_CpG.txt", 'w')
    o2 = open (f"{outfile_head}_CHG.txt", 'w')
    o3 = open (f"{outfile_head}_CHH.txt", 'w')

    margin = window//step

    for Chr in Chr_dic:
        Chr_lst = Chr_dic[Chr]

        for i in range(len(Chr_lst)-margin+1):
            
            count_lst = Chr_lst[i]
            index = i*step + window//2
            CpG_mC, CpG_C, CHG_mC, CHG_C, CHH_mC, CHH_C = count_lst
            CpG_total = CpG_mC + CpG_C
            CHG_total = CHG_mC + CHG_C
            CHH_total = CHH_mC + CHH_C

            if CpG_total == 0:
                CpG_perc = 0
            else:
                CpG_perc = round(Decimal(CpG_mC/(CpG_total) * 100),2)

            if CHG_total == 0:
                CHG_perc = 0
            else:
                CHG_perc = round(Decimal(CHG_mC/(CHG_total) * 100),2)
            
            if CHH_total == 0:
                CHH_perc = 0
            else:
                CHH_perc = round(Decimal(CHH_mC/(CHH_total) * 100),2)

            print(f"{Chr}\t{index}\t{CpG_mC}\t{CpG_C}\t{CpG_perc}", file = o1)
            print(f"{Chr}\t{index}\t{CHG_mC}\t{CHG_C}\t{CHG_perc}", file = o2)
            print(f"{Chr}\t{index}\t{CHH_mC}\t{CHH_C}\t{CHH_perc}", file = o3)
    
    o1.close()
    o2.close()
    o3.close()


def ProcessBP (bp):

    if bp >= 1000:
        text = f"{bp//1000}kb"
    else:
        text = f"{bp}bp"

    return text

    
def main():

    infile = sys.argv[1]
    Chr_file = sys.argv[2]
    outfile_head = sys.argv[3]
    cutoff = int(sys.argv[4])
    Window = int(sys.argv[5])
    Step = int(sys.argv[6])

    Window_text = ProcessBP (Window)
    Step_text   = ProcessBP (Step)

    outfile_filtered = f"{outfile_head}_filtered_{cutoff}.txt"
    outfile_comp_head = f"{outfile_head}_filtered_{cutoff}_compressed_{Window_text}_{Step_text}"
    
    Filter_CX_report (infile, outfile_filtered, cutoff)
    Chr_dic = Compress_CX_report (outfile_filtered, Chr_file, Window, Step)
    Print_Chr_dic (outfile_comp_head, Chr_dic, Window, Step)
    
    
main()
    
  
    
    

