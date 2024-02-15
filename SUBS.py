"""
Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.

01/02/2024
"""

# Imports
import os
import re

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def find_subs(s):
    s, t = s.split('\n')
    subs = re.finditer(f'(?=({t}))', s)
    locs = ' '.join([str(x.span()[0]+1) for x in subs])
    return locs

# Sample Test
sample_data = """GATATATGCATATACTT
ATAT"""
sample_output = find_subs(sample_data)
sample_answer = '2 4 10'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = '\n'.join([l.strip('\n') for l in f.readlines()])
    answer = find_subs(data)
    print(answer, '\n')
else:
    print('input file not found\n')
