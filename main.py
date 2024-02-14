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

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    aligner = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)

    hs_gg = [
        *aligner.align(hs_seq, gg_seq),
        "Gallus gallus",
    ]  # [0] is score, [1] is hs, [2] is gg, [3] is species name
    hs_mm = [*aligner.align(hs_seq, mm_seq), "Mus musculus"]
    hs_br = [*aligner.align(hs_seq, br_seq), "Balaeniceps rex"]
    hs_tt = [*aligner.align(hs_seq, tt_seq), "Tursiops truncatus"]

    alignments = [hs_gg, hs_mm, hs_br, hs_tt]
    alignments.sort(key=lambda x: x[0], reverse=True)

    print(
        "Alignments in order of similarity (first sequence is always Homo sapiens BRD2):"
    )
    for alignment in alignments:
        print(f"Homo sapiens vs. {alignment[3]}:\n{alignment[1]}\n{alignment[2]}\n")

    # print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("Scores after aligning all to Homo sapiens BRD2:")
    print(f"Gallus gallus: {hs_gg[0]}")
    print(f"Mus musculus: {hs_mm[0]}")
    print(f"Balaeniceps rex: {hs_br[0]}")
    print(f"Tursiops truncatus: {hs_tt[0]}")


if __name__ == "__main__":
    main()
