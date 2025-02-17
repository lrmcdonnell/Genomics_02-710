from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
import os

def run_msa(input_fasta, output_aln="aligned.fasta"):
    """
    Performs multiple sequence alignment using Clustal Omega.
    
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_aln (str): Path to store the aligned sequences (default: aligned.fasta).
    
    Returns:
        str: Path to the aligned output file.
    """
    clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta, outfile=output_aln, seqtype="DNA", verbose=True, auto=True)
    stdout, stderr = clustalomega_cline()
    return output_aln

def find_longest_consensus(aligned_fasta):
    """
    Finds the longest consensus sequence from a multiple sequence alignment.
    
    Parameters:
        aligned_fasta (str): Path to the aligned FASTA file.
    
    Returns:
        str: Longest consensus sequence.
    """
    alignment = AlignIO.read(aligned_fasta, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    
    # Get consensus sequence (threshold 0.7 conserves most common base at each position)
    consensus_seq = summary_align.dumb_consensus(threshold=0.7, ambiguous="N")
    
    # Convert to string and return longest continuous non-'N' sequence
    consensus_str = str(consensus_seq)
    longest_consensus = max(consensus_str.split("N"), key=len)
    
    return longest_consensus

if __name__ == "__main__":
    input_fasta = "sequences.fasta"  # Change this to your input FASTA file
    aligned_fasta = run_msa(input_fasta)
    longest_consensus = find_longest_consensus(aligned_fasta)
    
    print("Longest Consensus Sequence:\n", longest_consensus)
