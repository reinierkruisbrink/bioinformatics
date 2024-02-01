"""
Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u.
Given: A DNA string t having length at most 1000 nt.
Return: The transcribed RNA string of t.

01/02/2024
"""

# Imports
import os
import re

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def translate(s):
    return re.sub(r'[T]', 'U', s)

# Sample Test
sample_data = 'GATGGAACTTGACTACGTAAATT'
sample_output = translate(sample_data)
sample_answer = 'GAUGGAACUUGACUACGUAAAUU'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = f.readlines()[0].strip('\n')
    answer = translate(data)
    print(answer, '\n')
else:
    print('input file not found\n')
