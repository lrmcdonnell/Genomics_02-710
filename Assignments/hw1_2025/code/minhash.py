# minhash.py
# HW1, Computational Genomics, Spring 2024
# andrewid: lmcdonne

# WARNING: Do not change the file name; Autograder expects it.

import sys
import hashlib
import numpy as np

def ReadFASTA(filename):
    fp=open(filename, 'r')
    Sequences={}
    tmpname=""
    tmpseq=""
    for line in fp:
        if line[0]==">":
            if len(tmpseq)!=0:
                Sequences[tmpname]=tmpseq
            tmpname=line.strip().split()[0][1:]
            tmpseq=""
        else:
            tmpseq+=line.strip()
    Sequences[tmpname]=tmpseq
    fp.close()
    return Sequences

# You may define any helper functions for MinHash algorithm here

def get_kmers(sequence, k):
    """Extracts all k-mers from a given sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def hash_kmer(kmer):
    """Hashes a k-mer using SHA256 and returns an integer hash value."""
    return int(hashlib.sha256(kmer.encode()).hexdigest(), 16)


# Do not change this function signature
def minhash(seq1, seq2, k):
    """Approximates the Jaccard index between k-mers of seq1 and seq2 using MinHash
    Input:
        seq1, seq2: string
        k: int
    Returns: 1 items as so:
    the Jaccard index score (as a float)
    """
    num_hashes = 100
    kmers1 = get_kmers(seq1, k)
    kmers2 = get_kmers(seq2, k)
  
    hash_values1 = sorted([hash_kmer(kmer) for kmer in kmers1])[:num_hashes]
    hash_values2 = sorted([hash_kmer(kmer) for kmer in kmers2])[:num_hashes]

    minhash_set1 = set(hash_values1)
    minhash_set2 = set(hash_values2)
    
    intersection = len(minhash_set1 & minhash_set2)
    union = len(minhash_set1 | minhash_set2)
    
    return intersection / union if union > 0 else 0.0

if __name__=="__main__":
    Sequences=ReadFASTA(sys.argv[1])
    k = int(sys.argv[2])
    assert len(Sequences.keys())==2, "fasta file contains more than 2 sequences."
    seq1=Sequences[list(Sequences.keys())[0]]
    seq2=Sequences[list(Sequences.keys())[1]]

    score = minhash(seq1, seq2, k)

    print('Score: ', score)
