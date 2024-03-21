#import Bio
from collections import defaultdict

def composition(k=5, text='CAATCCAAC'):
    return [text[i:i+k] for i in range(len(text) - k + 1)]

def PathToGenome(text='ACCGA CCGAA CGAAG GAAGC AAGCT'):
    patterns = text.split()
    genome = patterns[0] + ''.join([pattern[-1] for pattern in patterns[1:]])
    return genome

def format_adj_list(adj_list):
    return '\n'.join([f'{k}: {" ".join(adj_list[k])}' for k in adj_list])

def OverlapGraph(text='ATGCG GCATG CATGC AGGCA GGCAT GGCAC'):
    kmers = text.split()
    adj_list = defaultdict(list)
    for kmer in kmers:
        for kmer2 in kmers:
            if kmer != kmer2 and kmer[1:] == kmer2[:-1]:
                adj_list[kmer].append(kmer2)
    return format_adj_list(adj_list)

def DeBruijnGraph(k=4, text='AAGATTCTCTAAGA'):
    kmers = composition(k, text)
    adj_list = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:k-1]
        sufix = kmer[1:]
        adj_list[prefix].append(sufix)
    return format_adj_list(adj_list)

def DeBruijnGraphFromKmers(kmers='GAGG CAGG GGGG GGGA CAGG AGGG GGAG'):
    adj_list = defaultdict(list)
    for kmer in kmers.split():
        prefix = kmer[:-1]
        sufix = kmer[1:]
        adj_list[prefix].append(sufix)
    return format_adj_list(adj_list)

def EulerianCycle(adj_list):
    pass