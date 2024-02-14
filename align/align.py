# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

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
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

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
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub
    
    # scoring function for matching base pairs (if matching -> 1, else -1)
    def score(self, b1: str, b2: str) -> int:
        s=-1
        if (b1==b2):
            s=1
        return s
    
    # helper function to complete alignment matrix
    def fill_M(self, i: int, j: int) -> float:
        M=self._align_matrix[i-1][j-1]+self.score(self._seqB[i-1], self._seqA[j-1])
        Ix=self._gapA_matrix[i-1][j-1]+self.score(self._seqB[i-1], self._seqA[j-1])
        Iy=self._gapB_matrix[i-1][j-1]+self.score(self._seqB[i-1], self._seqA[j-1])
        return max(M, Ix, Iy)

    # helper function to complete gapA matrix
    def fill_Ix(self, i: int, j: int) -> float:
        M=self._align_matrix[i-1][j]+self.gap_open+self.gap_extend
        Ix=self._gapA_matrix[i-1][j]+self.gap_extend
        return max(M, Ix)

    # helper function to complete gapB matrix
    def fill_Iy(self, i: int, j: int) -> float:
        M=self._align_matrix[i][j-1]+self.gap_open+self.gap_extend
        Iy=self._gapB_matrix[i][j-1]+self.gap_extend
        return max(M, Iy)
    
    # helper function that combines above three functions in a specific order 
    def fill_alignment_matrix(self, i: int, j: int) -> Tuple[float, float, float]:
        Mij=self.fill_M(i, j)
        Ixij=self.fill_Ix(i, j)
        Iyij=self.fill_Iy(i, j)
        return Mij, Ixij, Iyij

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
        self._seqA = seqA # seq1 (y)
        self._seqB = seqB # seq2 (x)
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        self._align_matrix=np.zeros((len(self._seqB)+1, len(self._seqA)+1))
        self._gapA_matrix=np.zeros((len(self._seqB)+1, len(self._seqA)+1))
        self._gapB_matrix=np.zeros((len(self._seqB)+1, len(self._seqA)+1))

        # initialize values

        # first row and col of align matrix are -inf

        # first row of align matrix
        for j in (range(len(self._seqA)+1)):
            self._align_matrix[0][j]=-np.inf

        # first column of align matrix
        for i in (range(len(self._seqB)+1)):
            self._align_matrix[i][0]=-np.inf

        # M(0,0)=0
        self._align_matrix[0][0]=0

        # first row of gapA
        for j in (range(len(self._seqA)+1)):
            self._gapA_matrix[0][j]=-np.inf

        # first column of gapA
        for i in (range(len(self._seqB)+1)):
            self._gapA_matrix[i][0]=self.gap_open+self.gap_extend*i

        # first column of gapB
        for i in (range(len(self._seqB)+1)):
            self._gapB_matrix[i][0]=-np.inf

        # first row of gapB
        for j in (range(len(self._seqA)+1)):
            self._gapB_matrix[0][j]=self.gap_open+self.gap_extend*j 

        # pass
        
        # TODO: Implement global alignment here
        
        # iterate through rows and columns and call helper functions to fill out the three matrices
        for i in range(1, len(self._seqB)+1):
            for j in range(1, len(self._seqA)+1):
                # print(i,j)
                self._align_matrix[i][j], self._gapA_matrix[i][j], self._gapB_matrix[i][j]=self.fill_alignment_matrix(i , j)
        
        
        # pass      		
        		    
        # TESTING BLOCK START #################################################
        # return self._align_matrix, self._gapA_matrix, self._gapB_matrix        
        # TESTING BLOCK END ###################################################
    
        return self._backtrace()


    # for the values across all three matrices at (i,j), get the max, and return an update vector that describes where to move
    def get_next_move(self, i, j):
        mat_values={"M":self._align_matrix[i][j], "Ix":self._gapA_matrix[i][j], "Iy":self._gapB_matrix[i][j]}
        max_key=max(mat_values, key=mat_values.get)
        # check if there are multiple max values
        max_key_list=[key for key,val in mat_values.items() if val==mat_values[max_key]]
        return max_key_list
    
    def get_alignment_score(self):
        prev_label_seq1=None
        prev_label_seq2=None
        alignment_score=0
        alignment_print=[]
        for s1, s2 in zip(self.seqA_align, self.seqB_align):
            # if (mtype=="c"):
            if (s1!="-" and s2!="-"):
                if (s1==s2):
                    # alignment_score+=match_val
                    alignment_score+=self.sub_dict[(s1,s2)]
                    alignment_print.append("|")
                else:
                    # alignment_score+=mismatch_val
                    alignment_score+=self.sub_dict[(s1,s2)]
                    alignment_print.append("·")
            else:
                if (s1=="-"):
                    if (prev_label_seq1!="gap"):
                        alignment_score+=self.gap_open+self.gap_extend
                        prev_label_seq1="gap"
                    else:
                        alignment_score+=self.gap_extend
                if (s2=="-"):
                    if (prev_label_seq2!="gap"):
                        alignment_score+=self.gap_open+self.gap_extend
                        prev_label_seq2="gap"
                    else:
                        alignment_score+=self.gap_extend
                alignment_print.append(" ")
            # else:


        return alignment_score, alignment_print, ''.join(self.seqA_align), ''.join(self.seqB_align)

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
        # starting positions are the current positions (bottom right corner) - pos variable is changed while start_pos is our record of the starting position
        start_pos=(len(self._seqB), len(self._seqA)) # -> curr_i=len(seq2), curr_j=len(seq1)
        pos=start_pos

        # alignment lists for sequences (they should be matching in length so we can draw 1 to 1 alignments)
        seqA_align=[]
        seqB_align=[]
        update_vectors={"M":(-1,-1), "Ix":(-1,0), "Iy":(0,-1)} # where M -> go up left, Ix -> go up, Iy -> left

        # method="highroad" 
        method="lowroad"

        # backtrace using the chosen method
        while (sum(pos)!=0):
            print(pos)
            max_key_list=self.get_next_move(pos[0], pos[1])

            # in the case of ties - for the highroad method, we prefer M; for the lowroad method, we prefer gaps
            if len(max_key_list)>1:
                if (method=="highroad"):
                    max_key=max_key_list[0]
                if (method=="lowroad"):
                    max_key=max_key_list[1]
            else:
                max_key=max_key_list[0]

            update_vector=update_vectors[max_key]

            pos=tuple(map(sum,zip(pos,update_vector)))
            if (max_key=="M"):
                seqB_align.append(self._seqB[pos[0]])
                seqA_align.append(self._seqA[pos[1]])
            elif (max_key=="Ix"):
                seqB_align.append(self._seqB[pos[0]])
                seqA_align.append("-")
            elif (max_key=="Iy"):
                seqA_align.append(self._seqA[pos[1]])
                seqB_align.append("-")
                
        # reverse lists
        seqA_align.reverse()
        seqB_align.reverse()

        self.seqA_align=seqA_align
        self.seqB_align=seqB_align

        self.alignment_score, alignment_print, self.seqA_align, self.seqB_align=self.get_alignment_score()

        # pass

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
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
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
