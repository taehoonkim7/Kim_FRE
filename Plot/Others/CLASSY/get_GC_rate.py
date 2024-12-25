import sys
from Bio import SeqIO # type: ignore


def get_GC_rate_from_sequence (seq):

    GC_count = 0
    for char in seq:
        if char in ("G","C"):
            GC_count += 1

    return GC_count/len(seq)


def get_GC_rate_by_window (seq, window):

    total_window = len(seq) - window + 1
    dic_GC = {}

    for i in range(total_window):
        GC_rate = get_GC_rate_from_sequence(seq[i:i+window])
        dic_GC[i] = GC_rate

    return dic_GC


def main():
    in_fasta = sys.argv[1]
    window = int(sys.argv[2])
    out_txt = sys.argv[3]

    with open (out_txt, "w") as o:
        print ("gene_id,position,GC_rate", file = o)

        for record in SeqIO.parse(in_fasta, "fasta"):
            dic_GC = get_GC_rate_by_window (record.seq, window)

            for i in dic_GC:
                print (f"{record.id},{i},{dic_GC[i]:.2f}", file = o)
        

if __name__ == "__main__":
    main()