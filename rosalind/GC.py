"""
Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. 
Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.

01/02/2024
"""

# Imports
import os
from utils import parse_FASTA

# Get the problem name
problem = __file__.split('\\')[-1][:-3].lower()
print(f"problem: {problem}\n")

# Solution
def gc_content(s):
    # Compute GC content
    GCs = {k:(v.count('G') + v.count('C')) / (len(v)+1e-6) * 100 for k, v in s.items()}
    # Return max GC content
    GC_content_argmax = max(GCs, key=GCs.get)
    GC_content_max = GCs[GC_content_argmax]
    return f'{GC_content_argmax}\n{GC_content_max:.6f}'

# Sample Test
sample_data = """>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
sample_output = gc_content(parse_FASTA(sample_data))
sample_answer = """Rosalind_0808
60.919540"""
assert sample_output == sample_answer
print('sample test passed\n')

# Output
data_file = f"data/rosalind_{problem}.txt"
if os.path.exists(data_file):
    with open(data_file) as f:
        data = ''.join(f.readlines())
    fasta_dict = parse_FASTA(data)
    answer = gc_content(fasta_dict)
    print(answer, '\n')
else:
    print('input file not found\n')
