# Importing Dependencies
import numpy as np
from typing import Tuple


# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """

    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(
            sub_matrix_file
        )  # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, "r") as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if "#" not in line.strip() and start is False:
                    residue_list = [
                        k for k in line.strip().upper().split(" ") if k != ""
                    ]
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(" ") if k != ""]
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(
                        line
                    ), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(
                            line[res_1]
                        )
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO

        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm

        Parameters:
                seqA: str
                        the first string to be aligned
                seqB: str
                        the second string to be aligned with seqA

        Returns:
                (alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
                        the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # Initialize matrix private attributes for use in alignment
        m, n = len(seqA), len(
            seqB
        )  # Seq A along rows (vertical), Seq B along columns (horizontal)

        # Main alignment matrix will store the max scores in each cell
        self._align_matrix = np.zeros(
            (m + 1, n + 1)
        )  # All should be size (m + 1, n + 1) to include gaps in 0th row and column
        self._align_matrix[0, 1:] = [-np.inf] * n  # Initialize proper scoring
        self._align_matrix[1:, 0] = [-np.inf] * m
        self._back = np.zeros((m + 1, n + 1), dtype=int)  # For main alignment backtrace

        # Initialize specific gap matrices (Ix and Iy) in the video
        self._gapA_matrix = np.full(
            (m + 1, n + 1), -np.inf
        )  # Initialize to negative infinity
        self._gapB_matrix = np.full(
            (m + 1, n + 1), -np.inf
        )  # Initialize to negative infinity

        # Initialize first row and column
        for i in range(0, m + 1):
            self._gapA_matrix[i, 0] = self.gap_open + i * self.gap_extend
        for j in range(0, n + 1):
            self._gapB_matrix[0, j] = self.gap_open + j * self.gap_extend

        # Fill matrices
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = self.sub_dict[(seqA[i - 1], seqB[j - 1])]

                # Diagonal move for match/mismatch
                diagonal_score = (
                    self._align_matrix[i - 1, j - 1] + match_score
                )  # Match score

                # Handling gaps in A (seqA gap extension or opening)
                if self._back[i, j - 1] == 1:
                    left_score_align = (
                        self._align_matrix[i, j - 1] + self.gap_extend
                    )  # Score for extending gap in A
                else:
                    left_score_align = (
                        self._align_matrix[i, j - 1] + self.gap_open + self.gap_extend
                    )
                left_score_A = (
                    self._gapA_matrix[i, j - 1] + self.gap_extend
                )  # If came from place where gap was highest
                self._gapA_matrix[i, j] = max(left_score_A, left_score_align)

                # Handling gaps in B (seqB gap extension or opening)
                if self._back[i - 1, j] == 2:
                    up_score_align = (
                        self._align_matrix[i - 1, j] + self.gap_extend
                    )  # Score for extending gap in B
                else:
                    up_score_align = (
                        self._align_matrix[i - 1, j] + self.gap_open + self.gap_extend
                    )
                up_score_B = (
                    self._gapB_matrix[i - 1, j] + self.gap_extend
                )  # If came from place where gap was highest
                self._gapB_matrix[i, j] = max(up_score_B, up_score_align)

                # Choose the best score
                max_score = max(
                    diagonal_score,
                    self._gapA_matrix[i, j],
                    self._gapB_matrix[i, j],
                )
                if max_score == diagonal_score:
                    self._back[i, j] = 0  # Diagonal move during backtrace
                elif max_score == self._gapA_matrix[i, j]:
                    self._back[i, j] = 1  # Left (gap in A)
                else:
                    self._back[i, j] = 2  # Up (gap in B)
                self._align_matrix[i, j] = max_score

        self.alignment_score = self._align_matrix[m, n]
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO

        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.

        Parameters:
                None

        Returns:
                (alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
                        the score and corresponding strings for the alignment of seqA and seqB
        """
        i, j = len(self._seqA), len(self._seqB)
        while i > 0 and j > 0:
            if self._back[i, j] == 0:  # Diagonal move
                self.seqA_align = self._seqA[i - 1] + self.seqA_align
                self.seqB_align = self._seqB[j - 1] + self.seqB_align
                i -= 1
                j -= 1
            elif self._back[i, j] == 1:  # Left move (gap in A)
                self.seqA_align = "-" + self.seqA_align
                self.seqB_align = self._seqB[j - 1] + self.seqB_align
                j -= 1  # Move left
            else:  # Up move (gap in B)
                self.seqA_align = self._seqA[i - 1] + self.seqA_align
                self.seqB_align = "-" + self.seqB_align
                i -= 1  # Move up

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(
        ".fa"
    ), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
