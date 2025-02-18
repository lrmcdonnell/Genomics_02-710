# normalization.py
# HW2, Computational Genomics, Spring 2025
# andrewid: lmcdonne

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np
import matplotlib.pyplot as plt

PER_MILLION = 1/1000000
PER_KILOBASE = 1/1000

# Do not change this function signature
def rpkm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    total_counts = np.sum(raw_counts, axis=0) 
    gene_lengths = gene_lengths[:, np.newaxis]
    
    # Calculate RPKM
    # (raw_counts * 1e9) because:  raw_counts/(gene_length/1000) divided by (total_counts/1e6)
    rpkm_matrix = (raw_counts * 1e9) / (gene_lengths * total_counts)
    return rpkm_matrix
   
# define any helper function here    


# Do not change this function signature
def size_factor(raw_counts):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    valid = np.all(raw_counts > 0, axis=1)
    gene_geom_mean = np.exp(np.mean(np.log(raw_counts[valid, :]), axis=1))

    ratios = raw_counts[valid, :] / gene_geom_mean[:, np.newaxis]
    size_factors = np.median(ratios, axis=0)
    
    # Normalize the counts by dividing each sample's counts by its size factor
    norm_counts = raw_counts / size_factors
    return norm_counts
    

if __name__=="__main__":
    raw_counts=np.loadtxt(sys.argv[1])
    gene_lengths=np.loadtxt(sys.argv[2])
    
    rpkm1=rpkm(raw_counts, gene_lengths)
    size_factor1=size_factor(raw_counts)

    # TODO: write plotting code here
    # For example, you might compare the distributions of normalized counts:
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.boxplot(rpkm1)
    plt.title("RPKM Normalized Counts")
    plt.xlabel("Samples")
    plt.ylabel("RPKM")
    
    plt.subplot(1, 2, 2)
    plt.boxplot(size_factor1)
    plt.title("Size Factor Normalized Counts")
    plt.xlabel("Samples")
    plt.ylabel("Normalized Counts")
    
    plt.tight_layout()
    plt.show()
