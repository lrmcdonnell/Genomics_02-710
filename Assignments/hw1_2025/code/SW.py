# SW.py
# HW1, Computational Genomics, Spring 2024
# andrewid: lmcdonne

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np

def ReadFASTA(filename):
    fp = open(filename, 'r')
    Sequences = {}
    tmpname = ""
    tmpseq = ""
    for line in fp:
        if line[0] == ">":
            if len(tmpseq) != 0:
                Sequences[tmpname] = tmpseq
            tmpname = line.strip().split()[0][1:]
            tmpseq = ""
        else:
            tmpseq += line.strip()
    Sequences[tmpname] = tmpseq
    fp.close()
    return Sequences

# Smith-Waterman Algorithm Implementation
def smith_waterman(seq1, seq2):
    match_score = 2
    mismatch_penalty = -2
    gap_penalty = -1
    
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring and traceback matrices
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    traceback = np.zeros((m + 1, n + 1), dtype=int)  # 0: stop, 1: diagonal, 2: up, 3: left
    
    max_score = 0
    max_pos = (0, 0)
    
    # Fill the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1, j] + gap_penalty
            insert = score_matrix[i, j - 1] + gap_penalty
            
            score_matrix[i, j] = max(0, match, delete, insert)
            
            if score_matrix[i, j] == match:
                traceback[i, j] = 1
            elif score_matrix[i, j] == delete:
                traceback[i, j] = 2
            elif score_matrix[i, j] == insert:
                traceback[i, j] = 3
            
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)
    
    # Traceback to get alignment
    aligned_seq1, aligned_seq2 = "", ""
    i, j = max_pos
    while traceback[i, j] != 0:
        if traceback[i, j] == 1:  # Diagonal (match/mismatch)
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i, j] == 2:  # Up (gap in seq2)
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif traceback[i, j] == 3:  # Left (gap in seq1)
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
    
    return max_score, aligned_seq1, aligned_seq2

if __name__ == "__main__":
    Sequences = ReadFASTA(sys.argv[1])
    assert len(Sequences.keys()) == 2, "FASTA file must contain exactly 2 sequences."
    
    seq1 = Sequences[list(Sequences.keys())[0]]
    seq2 = Sequences[list(Sequences.keys())[1]]
    
    score, align1, align2 = smith_waterman(seq1, seq2)
    # score, align1, align2 = smith_waterman('GGTTGACTA', 'TGTTACGG')
    print('Score:', score)
    print('Seq1:', align1)
    print('Seq2:', align2)
