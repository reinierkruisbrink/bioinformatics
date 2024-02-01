"""
Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t)
.

01/02/2024
"""

# Imports
import os

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def dist_hamm(s):
    s, t = s.split('\n')
    dist = sum([1 for a,b in zip(s,t) if a != b])
    return str(dist)

# Sample Test
sample_data = """GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT"""
sample_output = dist_hamm(sample_data)
sample_answer = '7'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = ''.join(f.readlines()).strip('\n')
    answer = dist_hamm(data)
    print(answer, '\n')
else:
    print('input file not found\n')
