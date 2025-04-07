import numpy as np

def calculate_motif_probability(sequence, pwm, nucleotides='ACGT'):
    """
    Calculate the probability of a given sequence according to a Position Weight Matrix (PWM).
    
    Parameters:
    -----------
    sequence : str
        The DNA sequence to calculate probability for
    pwm : dict or numpy.ndarray
        Position Weight Matrix with probabilities for each nucleotide at each position
    nucleotides : str
        String containing the nucleotides in the order they appear in the PWM
        
    Returns:
    --------
    float
        The probability of observing the given sequence based on the PWM
    """
    # Check if sequence length matches the PWM width
    if isinstance(pwm, dict):
        width = len(list(pwm.values())[0])
    else:
        width = pwm.shape[1]
    
    if len(sequence) != width:
        raise ValueError(f"Sequence length ({len(sequence)}) does not match PWM width ({width})")
    
    # Calculate the probability
    probability = 1.0
    
    for i, nucleotide in enumerate(sequence):
        # Convert nucleotide to index
        if nucleotide not in nucleotides:
            raise ValueError(f"Invalid nucleotide '{nucleotide}'. Must be one of {nucleotides}")
        
        nuc_idx = nucleotides.index(nucleotide)
        
        # Get probability from PWM
        if isinstance(pwm, dict):
            # If PWM is a dictionary with nucleotides as keys
            prob = pwm[nucleotide][i]
        else:
            # If PWM is a numpy array with rows corresponding to nucleotides
            prob = pwm[nuc_idx, i]
        
        # Multiply the probability
        probability *= prob
    
    return probability

# Define the PWM from the problem
pwm = np.array([
    [0.8, 0.1, 0.0, 0.9, 0.0, 0.3],  # A
    [0.0, 0.4, 0.05, 0.03, 0.1, 0.2], # C
    [0.2, 0.0, 0.95, 0.02, 0.1, 0.1], # G
    [0.0, 0.5, 0.0, 0.05, 0.8, 0.4]   # T
])

# Calculate the probability of sequence "ACCTTA"
sequence = "ACCTTA"
nucleotides = "ACGT"  # Order of nucleotides in the PWM

probability = calculate_motif_probability(sequence, pwm, nucleotides)
print(f"The probability of observing the motif {sequence} is: {probability}")
print(f"Step-by-step calculation:")
print("P(A at position 1) = 0.8")
print("P(C at position 2) = 0.4")
print("P(C at position 3) = 0.05")
print("P(T at position 4) = 0.05")
print("P(T at position 5) = 0.8")
print("P(A at position 6) = 0.3")
print(f"P(ACCTTA) = 0.8 × 0.4 × 0.05 × 0.05 × 0.8 × 0.3 = {probability}")