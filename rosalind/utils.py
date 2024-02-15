import re

def parse_FASTA(s):
    # s = s.split('>')[1:]
    ids = re.findall(r'>Rosalind_\d{1,4}', s)
    seqs = re.split(r'>Rosalind_\d{1,4}', s)
    fasta_dict = {}
    for i, seq in zip(ids, seqs[1:]):
        fasta_dict[i[1:]] = re.sub('\n', '', seq)
    return fasta_dict


if __name__ == "__main__":
    with open('data/rosalind_gc.txt') as f:
        data = ''.join(f.readlines())
        print(parse_FASTA(data))