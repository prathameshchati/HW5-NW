# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    # pass

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    # pass

    # create instance
    gap_open=-10
    gap_extend=-1
    nw=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open, gap_extend)

    # score dictionary    
    scores={}
    scores["Gallus_gallus"]=nw.align(hs_seq, gg_seq)[0]
    scores["Mus_musculus"]=nw.align(hs_seq, mm_seq)[0]
    scores["Balaeniceps_rex"]=nw.align(hs_seq, br_seq)[0]
    scores["tursiops_truncatu"]=nw.align(hs_seq, tt_seq)[0]

    # sort dictionary based on score
    scores={k: v for k, v in sorted(scores.items(), key=lambda s: s[1], reverse=True)}

    # print outputs
    print("BRD2 sequence similarity to Homo Sapiens:")
    for species, score in scores.items():
        print(f"{species}: {score}")

if __name__ == "__main__":
    main()
