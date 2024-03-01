"""
Given: A DNA string s of length at most 1000 nt.
Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.

01/02/2024
"""

# Imports
import os
from collections import Counter

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")


# Solution
def nt_counter(s):
    # output needs to be A C G T
    counter = Counter(s)
    return f"{counter['A']} {counter['C']} {counter['G']} {counter['T']}"

# Sample Test
sample_data = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
sample_output = nt_counter(sample_data)
sample_answer = '20 12 17 21'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = f.readlines()[0].strip('\n')
    answer = nt_counter(data)
    print(answer, '\n')
else:
    print('input file not found\n')