"""
Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.

01/02/2024
"""

# Imports
import os
from itertools import permutations

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def dom_pheno(s):
    k, m, n = map(int, s.split())
    pool = []
    pool += ['AA'] * k
    pool += ['Aa'] * m
    pool += ['aa'] * n
    pairs = [''.join(p) for p in permutations(pool, 2)]
    dom = 0
    for pair in pairs:
        if pair[1] == 'A' or pair[3] == 'A': # at least 1 homozygous dominant -> 100% chance of dominant
            dom += 1
        elif pair[0] == 'A' and pair[2] == 'A': # both heterozygous -> 75% chance of dominant
            dom += 0.75
        elif pair[0] == 'a' and pair[2] == 'a': # both homozygous recessive -> 0% chance of dominant
            continue
        else: # one heterozygous, one homozygous recessive -> 50% chance of dominant
           dom += 0.5
    dom_prob = f"{(dom / len(pairs)):.5f}"
    return dom_prob

# Sample Test
sample_data = '2 2 2'
sample_output = dom_pheno(sample_data)
sample_answer = '0.78333'
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = f.readlines()[0].strip('\n')
    answer = dom_pheno(data)
    print(answer, '\n')
else:
    print('input file not found\n')
