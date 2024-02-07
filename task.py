"""
task descriptions

dd/mm/2024
"""

# Imports
import os

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def task(s):
    return s

# Sample Test
sample_data = ''
sample_output = task(sample_data)
sample_answer = ''
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = f.readlines()[0].strip('\n')
    answer = task(data)
    print(answer, '\n')
else:
    print('input file not found\n')
