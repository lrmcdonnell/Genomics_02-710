def calculate_mec_score(reads, haplotype1, haplotype2):
    """
    Calculate the Minimum Error Correction (MEC) score for a set of reads and two haplotypes.
    
    Args:
        reads (list of str): List of reads, where '-' indicates no coverage at that position
        haplotype1 (str): First haplotype sequence
        haplotype2 (str): Second haplotype sequence
    
    Returns:
        int: The MEC score
    """
    total_mec = 0
    
    for read in reads:
        # Compare read to each haplotype
        errors_hap1 = 0
        errors_hap2 = 0
        
        for i, nucleotide in enumerate(read):
            if nucleotide == '-':
                continue  # Skip positions not covered by this read
            
            # Check if this position doesn't match haplotype 1
            if nucleotide != haplotype1[i]:
                errors_hap1 += 1
            
            # Check if this position doesn't match haplotype 2
            if nucleotide != haplotype2[i]:
                errors_hap2 += 1
        
        # Add the minimum number of errors to the total MEC score
        total_mec += min(errors_hap1, errors_hap2)
    
    return total_mec

# Input your reads here
reads = [
    "----------------1101011111-----------------------------",
    "----------10100100101000010011110----------------------",
    "-----001000101101101-----------------------------------",
    "0011000100010------------------------------------------",
    "------101110100100101----------------------------------",
    "------------------------------001010001010101----------",
    "------------------------------1101011101010------------",
    "----------------------------000010100010---------------",
    "----------------------00000011110101-------------------",
    "-----------------------------------------01010001010---",
    "---------------------------------------------1110101011",
    "-------------------------------------------010001110100",
    "-----------------------------------00010101010---------",
    "1100111010101------------------------------------------",
    "----------------001010000000---------------------------",
    "----------------------1111110010101--------------------",
    "-------------------------------------10101010111010----",
    "----------01011011010----------------------------------"
]

# Your proposed haplotypes
haplotype1 = "0011000100010110110101111111001010100010101010001110100"
haplotype2 = "1100111011101001001010000000110101011101010101110001011"

# Calculate MEC score
mec_score = calculate_mec_score(reads, haplotype1, haplotype2)
print(f"MEC Score: {mec_score}")

# If you want to see the comparison details for each read:
def show_read_details(reads, haplotype1, haplotype2):
    for i, read in enumerate(reads):
        errors_hap1 = 0
        errors_hap2 = 0
        
        for j, nucleotide in enumerate(read):
            if nucleotide == '-':
                continue
            
            if nucleotide != haplotype1[j]:
                errors_hap1 += 1
            
            if nucleotide != haplotype2[j]:
                errors_hap2 += 1
        
        min_errors = min(errors_hap1, errors_hap2)
        better_match = "Haplotype 1" if errors_hap1 <= errors_hap2 else "Haplotype 2"
        
        print(f"Read {i+1}: {errors_hap1} errors with Haplotype 1, {errors_hap2} errors with Haplotype 2")
        print(f"  Best match is {better_match} with {min_errors} errors")

# Uncomment the line below if you want to see details for each read
show_read_details(reads, haplotype1, haplotype2)