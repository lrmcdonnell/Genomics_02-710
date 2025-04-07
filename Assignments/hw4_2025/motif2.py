import numpy as np

def find_consensus_sequence(pwm, nucleotides='ACGT'):
    """
    Find the consensus sequence from a Position Weight Matrix (PWM).
    The consensus sequence contains the most probable nucleotide at each position.
    
    Parameters:
    -----------
    pwm : numpy.ndarray
        Position Weight Matrix with probabilities for each nucleotide at each position
    nucleotides : str
        String containing the nucleotides in the order they appear in the PWM
        
    Returns:
    --------
    str
        The consensus sequence
    float
        The probability of the consensus sequence
    """
    # Find the index of the maximum probability for each position
    max_indices = np.argmax(pwm, axis=0)
    
    # Convert indices to nucleotides
    consensus = ''.join(nucleotides[idx] for idx in max_indices)
    
    # Calculate the probability of the consensus sequence
    probability = 1.0
    for i, idx in enumerate(max_indices):
        probability *= pwm[idx, i]
    
    return consensus, probability

# Define the PWM from the problem
pwm = np.array([
    [0.8, 0.1, 0.0, 0.9, 0.0, 0.3],  # A
    [0.0, 0.4, 0.05, 0.03, 0.1, 0.2], # C
    [0.2, 0.0, 0.95, 0.02, 0.1, 0.1], # G
    [0.0, 0.5, 0.0, 0.05, 0.8, 0.4]   # T
])

nucleotides = "ACGT"  # Order of nucleotides in the PWM

# Find the consensus sequence
consensus, probability = find_consensus_sequence(pwm, nucleotides)

print(f"The consensus sequence is: {consensus}")
print(f"The probability of this sequence is: {probability}")

# Print the step-by-step identification of the consensus
print("\nStep-by-step identification of consensus:")
for i in range(pwm.shape[1]):
    position_probs = pwm[:, i]
    max_idx = np.argmax(position_probs)
    max_nucleotide = nucleotides[max_idx]
    max_prob = position_probs[max_idx]
    
    print(f"Position {i+1}: {nucleotides[0]}={position_probs[0]:.2f}, {nucleotides[1]}={position_probs[1]:.2f}, "
          f"{nucleotides[2]}={position_probs[2]:.2f}, {nucleotides[3]}={position_probs[3]:.2f} → "
          f"Max is {max_nucleotide} with probability {max_prob:.2f}")