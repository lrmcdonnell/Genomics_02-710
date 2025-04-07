import numpy as np
import matplotlib.pyplot as plt

def calculate_entropy(position_probs):
    """
    Calculate the entropy of a position in the PWM.
    Entropy H = -Σ[p_i * log2(p_i)]
    
    Parameters:
    -----------
    position_probs : numpy.ndarray
        Array of probabilities for each nucleotide at a specific position
        
    Returns:
    --------
    float
        The entropy value
    """
    # Replace 0 probabilities with a small value to avoid log(0)
    position_probs = np.where(position_probs > 0, position_probs, 1e-10)
    
    # Calculate entropy
    entropy = -np.sum(position_probs * np.log2(position_probs))
    return entropy

def calculate_information_content(entropy):
    """
    Calculate the information content of a position.
    R = log2(4) - H
    Note: This simplified version doesn't include the small-sample correction (e_n)
    
    Parameters:
    -----------
    entropy : float
        Entropy value
        
    Returns:
    --------
    float
        The information content
    """
    return np.log2(4) - entropy  # log2(4) = 2 bits

# Define the PWM from the problem
pwm = np.array([
    [0.8, 0.1, 0.0, 0.9, 0.0, 0.3],  # A
    [0.0, 0.4, 0.05, 0.03, 0.1, 0.2], # C
    [0.2, 0.0, 0.95, 0.02, 0.1, 0.1], # G
    [0.0, 0.5, 0.0, 0.05, 0.8, 0.4]   # T
])

# Calculate entropy for each position
entropies = []
info_contents = []

print("Position-by-position entropy calculation:")
print("----------------------------------------")
for i in range(pwm.shape[1]):
    position_probs = pwm[:, i]
    
    # Calculate and store entropy
    entropy = calculate_entropy(position_probs)
    entropies.append(entropy)
    
    # Calculate and store information content
    info_content = calculate_information_content(entropy)
    info_contents.append(info_content)
    
    print(f"Position {i+1}:")
    print(f"  Probabilities: A={position_probs[0]:.2f}, C={position_probs[1]:.2f}, G={position_probs[2]:.2f}, T={position_probs[3]:.2f}")
    print(f"  Entropy: {entropy:.4f} bits")
    print(f"  Information Content: {info_content:.4f} bits")
    print()

# Find position with greatest entropy
max_entropy_pos = np.argmax(entropies) + 1  # +1 for 1-based indexing
min_entropy_pos = np.argmin(entropies) + 1  # +1 for 1-based indexing

print(f"Position with greatest entropy: Position {max_entropy_pos} (Entropy = {max(entropies):.4f} bits)")
print(f"Position with lowest entropy: Position {min_entropy_pos} (Entropy = {min(entropies):.4f} bits)")

# Visualize the entropies and information content
plt.figure(figsize=(10, 6))
positions = np.arange(1, pwm.shape[1] + 1)

plt.subplot(2, 1, 1)
plt.bar(positions, entropies, color='blue', alpha=0.7)
plt.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='Maximum entropy (2 bits)')
plt.xlabel('Position')
plt.ylabel('Entropy (bits)')
plt.title('Entropy by Position')
plt.xticks(positions)
plt.ylim(0, 2.1)
plt.legend()

plt.subplot(2, 1, 2)
plt.bar(positions, info_contents, color='green', alpha=0.7)
plt.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='Maximum information content (2 bits)')
plt.xlabel('Position')
plt.ylabel('Information Content (bits)')
plt.title('Information Content by Position')
plt.xticks(positions)
plt.ylim(0, 2.1)
plt.legend()

plt.tight_layout()
plt.savefig('entropy_info_content.png', dpi=300, bbox_inches='tight')

# Print relationship between entropy and conservation
print("\nRelationship between entropy and conservation:")
print("----------------------------------------------")
print("1. Low entropy indicates high conservation: When a position strongly prefers one nucleotide over others.")
print("2. High entropy indicates low conservation: When a position has more uniform distribution of nucleotides.")
print("3. Maximum entropy (2 bits) occurs when all four nucleotides are equally likely (0.25 each).")
print("4. Minimum entropy (0 bits) occurs when one nucleotide has probability 1 and others have 0.")
print("5. Information content is inversely related to entropy: High information content = Low entropy = High conservation.")