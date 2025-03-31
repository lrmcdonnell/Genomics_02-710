#!/usr/bin/env python3

import os
from collections import Counter, defaultdict


def read_fasta(fasta_file):
    """
    Read a FASTA file and return a dictionary of sequences.
    
    Args:
        fasta_file (str): Path to the FASTA file
        
    Returns:
        dict: Dictionary with header as key and sequence as value
    """
    sequences = {}
    current_header = None
    current_sequence = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = current_sequence
                current_header = line[1:]  # Remove the '>' character
                current_sequence = ""
            else:
                current_sequence += line
                
    if current_header and current_sequence:
        sequences[current_header] = current_sequence
        
    return sequences


def extract_kmers(sequence, k):
    """
    Extract all k-mers from a sequence.
    
    Args:
        sequence (str): DNA sequence
        k (int): k-mer length
        
    Returns:
        set: Set of k-mers in the sequence
    """
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if all(base.upper() in 'ACGT' for base in kmer):
            kmers.add(kmer)
    return kmers


def count_kmers(sequence, k):
    """
    Count k-mers in a sequence.
    
    Args:
        sequence (str): DNA sequence
        k (int): k-mer length
        
    Returns:
        Counter: Count of each k-mer
    """
    kmers = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if all(base.upper() in 'ACGT' for base in kmer):
            kmers[kmer] += 1
    return kmers


def find_discriminative_kmers(data_folder, k=4):
    """
    Find discriminative k-mers for each genome.
    
    Args:
        data_folder (str): Path to the folder containing genome FASTA files
        k (int): k-mer length
        
    Returns:
        dict: Dictionary with genome filenames as keys and their discriminative k-mers as values
    """
    fasta_files = [f for f in os.listdir(data_folder) 
                  if (f.endswith('.fasta') or f.endswith('.fa')) and f != 'reads.fa']
    
    print(f"Found {len(fasta_files)} genome files: {fasta_files}")
    
    # Extract k-mers
    genome_kmers = {}
    kmer_occurrences = defaultdict(set) 
    
    for fasta_file in fasta_files:
        print(f"Extracting k-mers from genome: {fasta_file}")
        
        # Read FASTA 
        file_path = os.path.join(data_folder, fasta_file)
        sequences = read_fasta(file_path)
        
        sequence = ''
        for header, seq in sequences.items():
            sequence += seq
        
        print(f"Sequence length: {len(sequence)} bp")
        
        kmer_set = extract_kmers(sequence, k)
        kmer_counts = count_kmers(sequence, k)
        
        print(f"Found {len(kmer_set)} unique {k}-mers in {fasta_file}")
        
        genome_kmers[fasta_file] = {'kmers': kmer_set, 'counts': kmer_counts}
        
        for kmer in kmer_set:
            kmer_occurrences[kmer].add(fasta_file)
    
   
    kmer_genome_counts = Counter([len(genomes) for genomes in kmer_occurrences.values()])
    for count, num_kmers in sorted(kmer_genome_counts.items()):
        print(f"K-mers appearing in exactly {count} genomes: {num_kmers}")
    
   
    discriminative_kmers = {}
    
    for target_genome in fasta_files:
        print(f"\nFinding discriminative k-mers for genome: {target_genome}")
        
        target_kmers = genome_kmers[target_genome]['kmers']
        target_kmer_counts = genome_kmers[target_genome]['counts']
        
        unique_kmers = set()
        for kmer in target_kmers:
            if len(kmer_occurrences[kmer]) == 1 and target_genome in kmer_occurrences[kmer]:
                unique_kmers.add(kmer)
        
        # counts for discriminative k-mers
        discriminative_kmer_counts = {kmer: target_kmer_counts[kmer] for kmer in unique_kmers}
        
        top_discriminative_kmers = sorted(discriminative_kmer_counts.items(), 
                                         key=lambda x: x[1], reverse=True)
        
        discriminative_kmers[target_genome] = top_discriminative_kmers
        
        print(f"Found {len(unique_kmers)} discriminative k-mers")
        if top_discriminative_kmers:
            print(f"Top 5 most common discriminative {k}-mers:")
            for i, (kmer, count) in enumerate(top_discriminative_kmers[:5], 1):
                print(f"{i}. {kmer}: {count}")
    
    return discriminative_kmers


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
    
    print(f"Using {len(kmer_to_genome)} discriminative k-mers for classification")
    
    reads = read_fasta(reads_file)
    
    classifications = Counter()
    unclassified = 0
    total_reads = len(reads)
    
    print(f"Classifying {total_reads} reads using discriminative k-mers...")
    
    for i, (header, read) in enumerate(reads.items()):
        if (i+1) % 100 == 0:
            print(f"Processed {i+1}/{total_reads} reads")
            
        read_kmers = extract_kmers(read, k)
        
        genome_matches = Counter()
        
        for kmer in read_kmers:
            if kmer in kmer_to_genome:
                genome = kmer_to_genome[kmer]
                genome_matches[genome] += 1
        
        # Classify based on the genome with the most matching discriminative k-mers
        if genome_matches:
            best_match = genome_matches.most_common(1)[0][0]
            classifications[best_match] += 1
        else:
            unclassified += 1
    
    frequencies = {genome: count/total_reads for genome, count in classifications.items()}
    
    if unclassified > 0:
        classifications['unclassified'] = unclassified
        frequencies['unclassified'] = unclassified / total_reads
    
    return classifications, frequencies


def main():
    k = 4
    data_folder = "data"
    
    print("Finding discriminative k-mers...")
    discriminative_kmers = find_discriminative_kmers(data_folder, k)
    
    reads_file = os.path.join(data_folder, "reads.fa")
    classifications, frequencies = classify_reads_with_discriminative_kmers(
        reads_file, discriminative_kmers, k)
    
    print("\nClassification Results:")
    print("=======================")
    
    for genome, count in sorted(classifications.items()):
        if genome != 'unclassified':
            print(f"Genome: {genome}")
            print(f"Count: {count}")
            print(f"Relative Frequency: {frequencies[genome]:.4f} ({frequencies[genome]*100:.2f}%)")
            print()
    
    if 'unclassified' in classifications:
        print(f"Unclassified reads: {classifications['unclassified']}")
        print(f"Unclassified frequency: {frequencies['unclassified']:.4f} ({frequencies['unclassified']*100:.2f}%)")


if __name__ == "__main__":
    main()