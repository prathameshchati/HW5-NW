# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np
import Bio
from Bio import pairwise2


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

    # compute gapA, gapB and alignment matrix and see if the outputs of function match
    gapA_matrix=np.array([[-10., -np.inf, -np.inf, -np.inf, -np.inf],
                          [-11., -22., -23., -24., -25.],
                          [-12.,  -6., -23., -23., -25.],
                          [-13.,  -7.,  -7., -12., -17.]])
    gapB_matrix=np.array([[-10., -11., -12., -13., -14.],
                          [-np.inf, -22.,  -6.,  -7.,  -8.],
                          [-np.inf, -23., -22.,  -7.,  -8.],
                          [-np.inf, -24., -24., -19.,  -6.]])
    align_matrix=np.array([[  0., -11., -12., -13., -14.],
                           [-11.,   5., -12., -12., -14.],
                           [-12., -11.,   4.,  -1.,  -6.],
                           [-13., -13.,  -8.,   5.,   4.]])
    
    # initialize nw with gap open and extend penalties and blosum62 matrix
    gap_open=-10
    gap_extend=-1
    nw=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open, gap_extend)

    # align
    alignment_score, seqA_align, seqB_align=nw.align(seq1, seq2)

    # use pairwise2 module from Biopython package to also align and check if correct solution is given # https://biopython.org/docs/1.75/api/Bio.Align.html
    alignments=pairwise2.align.globalds(seq1, seq2, nw.sub_dict, gap_open+gap_extend, gap_extend)
    pw2_score=alignments[0].score
    pw2_seqA=alignments[0].seqA
    pw2_seqB=alignments[0].seqB

    # check that three matrices are being correctly outputted
    assert np.array_equal(nw._gapA_matrix, gapA_matrix)
    assert np.array_equal(nw._gapB_matrix, gapB_matrix)
    assert np.array_equal(nw._align_matrix, align_matrix)

    # assert alingment score is correct and at the bottom right value of the alignment matrix and that it matches the pairwise2 outputs
    assert alignment_score==pw2_score
    assert nw._align_matrix[len(seq2)][len(seq1)]==pw2_score

    # assert sequence alignments are the same
    assert seqA_align==pw2_seqA
    assert seqB_align==pw2_seqB

    # pass

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    # initialize nw with gap open and extend penalties and blosum62 matrix
    gap_open=-10
    gap_extend=-1
    nw=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open, gap_extend)

    # align
    alignment_score, seqA_align, seqB_align=nw.align(seq3, seq4)

    # use pairwise2 module from Biopython package to also align and check if correct solution is given # https://biopython.org/docs/1.75/api/Bio.Align.html
    alignments=pairwise2.align.globalds(seq3, seq4, nw.sub_dict, gap_open+gap_extend, gap_extend)
    pw2_score=alignments[0].score
    pw2_seqA=alignments[0].seqA
    pw2_seqB=alignments[0].seqB

    # if backtracing is done correctly, the alignments of the sequence should match the exact outputs of pairwise2 (with the specific gaps and matches)
    # this way we can check what we tested above
    
    # assert alingment score is correct and at the bottom right value of the alignment matrix and that it matches the pairwise2 outputs
    assert alignment_score==pw2_score # score is 17 
    assert nw._align_matrix[len(seq4)][len(seq3)]==pw2_score # score is 17 

    # assert sequence alignments are the same
    assert seqA_align==pw2_seqA
    assert seqB_align==pw2_seqB

    # pass


# check that the brd2 sequences are yielding the correct alignment scores when compared to the human sequence
def test_large_sequences():

    # read fasta sequences
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # create instance
    gap_open=-10
    gap_extend=-1
    nw=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open, gap_extend)

    # score
    hs_gg=nw.align(hs_seq, gg_seq)
    hs_mm=nw.align(hs_seq, mm_seq)
    hs_br=nw.align(hs_seq, br_seq)
    hs_tt=nw.align(hs_seq, tt_seq)

    # run pairwise2
    hs_gg_pw2=pairwise2.align.globalds(hs_seq, gg_seq, nw.sub_dict, gap_open+gap_extend, gap_extend)
    hs_mm_pw2=pairwise2.align.globalds(hs_seq, mm_seq, nw.sub_dict, gap_open+gap_extend, gap_extend)
    hs_br_pw2=pairwise2.align.globalds(hs_seq, br_seq, nw.sub_dict, gap_open+gap_extend, gap_extend)
    hs_tt_pw2=pairwise2.align.globalds(hs_seq, tt_seq, nw.sub_dict, gap_open+gap_extend, gap_extend)

    # we can check scores to make sure they are equal - however, since pw2 produces numerous alignments (longer sequences being used), the exact alignment is not being checked
    assert hs_gg[0]==hs_gg_pw2[0].score
    assert hs_mm[0]==hs_mm_pw2[0].score
    assert hs_br[0]==hs_br_pw2[0].score
    assert hs_tt[0]==hs_tt_pw2[0].score

    # pass




