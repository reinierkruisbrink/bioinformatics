"""
Given: Positive integers n≤40 and k≤5.

Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).

01/02/2024
"""

# Imports
import os

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def reproduction(s):
    n, k = map(int, s.split())
    return str(fib(n, k))

def fib(n, k):
    if n == 1 or n == 2: # base case
        return 1
    else: # recursive case (prev fib + litter size * prev prev fib)
        return fib(n-1, k) + k*fib(n-2, k)

# Sample Test
sample_data = '5 3'
sample_output = reproduction(sample_data)
sample_answer = '19'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = f.readlines()[0].strip('\n')
    answer = reproduction(data)
    print(answer, '\n')
else:
    print('input file not found\n')
