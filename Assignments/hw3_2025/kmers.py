#!/usr/bin/env python3

import os
from collections import Counter, defaultdict
from discriminative_kmers import read_fasta, find_discriminative_kmers



def extract_read_kmers(read, k):
    """
    Extract all k-mers from a read.
    
    Args:
        read (str): DNA read sequence
        k (int): k-mer length
        
    Returns:
        set: Set of k-mers in the read
    """
    kmers = set()
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        if all(base in 'ACGT' for base in kmer):
            kmers.add(kmer)
    return kmers


def build_discriminative_kmer_index(discriminative_kmers):
    """
    Build an index mapping discriminative k-mers to their genome of origin.
    
    Args:
        discriminative_kmers (dict): Dictionary with genome filenames as keys 
                                    and their discriminative k-mers as values
        
    Returns:
        dict: Dictionary mapping each discriminative k-mer to its genome
    """
    kmer_to_genome = {}
    
    for genome, kmers in discriminative_kmers.items():
        for kmer, count in kmers:
            kmer_to_genome[kmer] = genome
    
    return kmer_to_genome


def classify_reads_with_discriminative_kmers(reads_file, discriminative_kmers, k=10):
    """
    Classify reads using discriminative k-mers.
    
    Args:
        reads_file (str): Path to the reads FASTA file
        discriminative_kmers (dict): Dictionary with genome filenames as keys 
                                    and their discriminative k-mers as values
        k (int): k-mer length
        
    Returns:
        tuple: (classifications, frequencies) - counts and relative frequencies of reads per genome
    """
    kmer_to_genome = build_discriminative_kmer_index(discriminative_kmers)
    reads = read_fasta(reads_file)
    
    # Classify each read
    classifications = Counter()
    unclassified = 0
    total_reads = len(reads)
    
    print(f"Classifying {total_reads} reads using discriminative k-mers...")
    
    for i, (header, read) in enumerate(reads.items()):
        if (i+1) % 100 == 0:
            print(f"Processed {i+1}/{total_reads} reads")
            
        # Extract k-mers
        read_kmers = extract_read_kmers(read, k)
        
        genome_matches = Counter()
        
        for kmer in read_kmers:
            if kmer in kmer_to_genome:
                genome = kmer_to_genome[kmer]
                genome_matches[genome] += 1
        
        if genome_matches:
            best_match = genome_matches.most_common(1)[0][0]
            classifications[best_match] += 1
        else:
            unclassified += 1
    
    total_classified = sum(classifications.values())
    frequencies = {genome: count/total_reads for genome, count in classifications.items()}

    if unclassified > 0:
        classifications['unclassified'] = unclassified
        frequencies['unclassified'] = unclassified / total_reads
    
    return classifications, frequencies


def main():
    k = 10
    data_folder = "data"
    
    discriminative_kmers = find_discriminative_kmers(data_folder, k)

    reads_file = os.path.join(data_folder, "reads.fa")
    classifications, frequencies = classify_reads_with_discriminative_kmers(
        reads_file, discriminative_kmers, k)
    
    print("\nClassification Results using Discriminative K-mers:")
    print("==================================================")
    
    for genome, count in classifications.items():
        print(f"Genome: {genome}")
        print(f"Count: {count}")
        print(f"Relative Frequency: {frequencies[genome]:.4f} ({frequencies[genome]*100:.2f}%)")
        print()

if __name__ == "__main__":
    main()