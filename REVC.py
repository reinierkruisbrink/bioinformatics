"""
Given: A DNA string s of length at most 1000 bp.
Return: The reverse complement sc of s.

01/02/2024
"""

# Imports
import os
import re

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def complement(s):
    return re.sub(r'[ATCG]', lambda x: {'A':'T', 'T':'A', 'C':'G', 'G':'C'}[x.group(0)], s[::-1])

# Sample Test
sample_data = 'AAAACCCGGT'
sample_output = complement(sample_data)
sample_answer = 'ACCGGGTTTT'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = f.readlines()[0].strip('\n')
    answer = complement(data)
    print(answer, '\n')
else:
    print('input file not found\n')
