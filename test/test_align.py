# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np


def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    aligner = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    alignment = aligner.align(seq1, seq2)

    assert alignment[1] == "MYQR"
    assert alignment[2] == "M-QR"
    assert aligner._align_matrix == [
        [0, -np.inf, -np.inf, -np.inf, -np.inf],
        [-np.inf, 5, -12, -12, -14],
        [-np.inf, -11, 4, 0, -5],
        [-np.inf, -13, -7, 5, 5],
    ]
    assert aligner._gapA_matrix == [
        [-10, -np.inf, -np.inf, -np.inf, -np.inf],
        [-11, -np.inf, -np.inf, -np.inf, -np.inf],
        [-12, -5, -22, -22, -24],
        [-13, -6, -6, -10, -15],
    ]
    assert aligner._gapB_matrix == [
        [-10, -11, -12, -13, -14],
        [-np.inf, -np.inf, -5, -6, -7],
        [-np.inf, -np.inf, -21, -6, -10],
        [-np.inf, -np.inf, -23, -17, -5],
    ]


def test_nw_backtrace():
    """
    Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    aligner = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    alignment = aligner.align(seq3, seq4)

    assert alignment[0] == 17
    assert alignment[1] == "MAVHQLIRRP"
    assert alignment[2] == "M---QLIRHP"
