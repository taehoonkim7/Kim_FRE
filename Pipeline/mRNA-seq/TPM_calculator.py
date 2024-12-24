import sys

def Make_dic_transcript_to_gene (f_trans_gene):
    
    f = open (f_trans_gene, 'r')
    ls = f.readlines()
    f.close()

    dic = {}
    for l in ls:
        es = l.strip().split('\t')
        gene = es[0]
        transcript = es[1]

        dic[transcript] = gene

    return dic


def Process_counts_table_to_dic (f_count):
    
    f = open (f_count, 'r')
    ls = f.readlines()
    f.close()

    new_ls = []
    for l in ls:
        new_ls.append (l.strip())
    
    l_samples = new_ls[0]
    counts_matrix = new_ls[1:]
    
    lst_samples = l_samples.strip().split('\t')[1:]
    dic_counts = {}
    lst_genes = []
    
    for sample in lst_samples:
        dic_counts[sample] = []
    
    for i in range(len(counts_matrix)):
        
        l = counts_matrix[i]
        elements = l.strip().split("\t")
        
        GeneID = elements[0]
        counts = elements[1:]
        
        lst_genes.append (GeneID)
        
        for j in range (len(counts)):
            dic_counts[lst_samples[j]].append (counts[j])
    
    return lst_samples, lst_genes, dic_counts
    

def Process_gff_to_dic_genelen (f_gff, dic_trans_gene):

    f = open (f_gff, 'r')
    lst_gff = f.readlines()
    f.close()
    lst_exon_gff = [line for line in lst_gff if "exon" in line]
    
    lst_transcript = list(dic_trans_gene.keys())
    
    dic_Transcript_gff = {}
    
    for line in lst_exon_gff:
        es = line.strip().split('\t')
        attr = es[8]
        ID = attr.split(";")[1].split("Parent=")[1]
        
        if ID in lst_transcript:
            if ID in dic_Transcript_gff:
                dic_Transcript_gff[ID].append (line)
            else:
                dic_Transcript_gff[ID] = [line]
    
    dic_GeneLen = {}
    
    for Transcript in dic_Transcript_gff:
        
        lst = dic_Transcript_gff [Transcript]
        Gene = dic_trans_gene [Transcript]
        
        length = 0
        for line in lst:
            elements = line.strip().split('\t')
            start = int(elements[3])
            end = int(elements[4])
            range = end - start + 1
            
            length += range
           
        dic_GeneLen [Gene] = length/1000
    
    return dic_GeneLen
    

def Calculate_TPM (lst_genes, dic_counts, dic_GeneLen):
    
    dic_TPM = {}
    
    for sample in dic_counts:
        
        lst_temp_count = dic_counts[sample]
        lst_RPK = []
        
        for i in range(len(lst_temp_count)):
            count = int(lst_temp_count[i])
            Gene = lst_genes[i]
            GeneLen = dic_GeneLen[Gene]
            RPK = count/GeneLen
            
            lst_RPK.append (RPK)
            
        sum_RPK = sum(lst_RPK)    
        
        lst_TPM = [str(RPK * (1000000/sum_RPK)) for RPK in lst_RPK]
        dic_TPM [sample] = lst_TPM
        
    return dic_TPM
    

def Format_TPM_matrix (lst_samples, lst_genes, dic_TPM):

    lst_lst_element = []
    
    for gene in lst_genes:
        lst_lst_element.append ([gene])
    
    for sample in lst_samples:
        lst_TPM = dic_TPM [sample]
        
        for i in range(len(lst_TPM)):
            lst_lst_element[i].append (lst_TPM[i])
            
    header = [["GeneID"] + lst_samples]
    lst_lst_element = header + lst_lst_element

    lst_lines =[]
    
    for lst in lst_lst_element:
        lst_lines.append ('\t'.join (lst))
        
    return lst_lines
    
    
def Print_output (outfilename, lst_lines):

    outfile = open (outfilename, 'w')
    for line in lst_lines:
        print (line, file = outfile)
    outfile.close()
    
    
def main():
    
    f_count = sys.argv[1] #count matrix (.tsv)
    f_out = sys.argv[2] #outfilename: TPM matrix (.tsv)
    f_gff3 = sys.argv[3] #gff3 file (.gff3)
    f_trans_gene = sys.argv[4] #list of primary transcripts and gene IDs (.txt) - tab delimited
    
    Lst_samples, Lst_genes, Dic_counts_per_sample = Process_counts_table_to_dic (f_count)
    
    dic_trans_gene = Make_dic_transcript_to_gene (f_trans_gene)
    Dic_GeneLen = Process_gff_to_dic_genelen (f_gff3, dic_trans_gene)
    
    Dic_TPM = Calculate_TPM (Lst_genes, Dic_counts_per_sample, Dic_GeneLen)
    TPM_matrix = Format_TPM_matrix (Lst_samples, Lst_genes, Dic_TPM)
    Print_output (f_out, TPM_matrix)
    
    
if __name__ == "__main__":
    main()