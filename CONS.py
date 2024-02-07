"""
Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)

07/02/2024
"""

# Imports
import os
import numpy as np
from utils import parse_FASTA

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def profile(s):
    # Only use the DNA sequences
    s = list(s.values())
    # Initialize profile matrix
    profiles = np.zeros((4, len(s[0]))).astype(int)
    # Count bases per sequence
    for seq in s:
        for i, base in enumerate(seq):
            if base == 'A':
                profiles[0, i] += 1
            elif base == 'C':
                profiles[1, i] += 1
            elif base == 'G':
                profiles[2, i] += 1
            elif base == 'T':
                profiles[3, i] += 1
    # Output formatting
    consensus = ''.join([{0:'A', 1:'C', 2:'G', 3:'T'}[arg] for arg in np.argmax(profiles, axis=0)])
    output = consensus + '\n'
    output += 'A: ' + ' '.join(profiles[0].astype(str)) + '\n'
    output += 'C: ' + ' '.join(profiles[1].astype(str)) + '\n'
    output += 'G: ' + ' '.join(profiles[2].astype(str)) + '\n'
    output += 'T: ' + ' '.join(profiles[3].astype(str))
    return output

# Sample Test
sample_data = """>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT"""
sample_output = profile(parse_FASTA(sample_data))
sample_answer = """ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6"""
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = ''.join(f.readlines())
    fasta_dict = parse_FASTA(data)
    answer = profile(fasta_dict)
    print(answer, '\n')
    # Use the output file to submit to rosalind site
    with open('output_file.txt', 'w') as f:
        f.write(answer)
else:
    print('input file not found\n')
