import random
import itertools
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
    return adj_list

def DeBruijnGraph(k=4, text='AAGATTCTCTAAGA'):
    kmers = composition(k, text)
    adj_list = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:k-1]
        sufix = kmer[1:]
        adj_list[prefix].append(sufix)
    return adj_list

def DeBruijnGraphFromKmers(kmers='GAGG CAGG GGGG GGGA CAGG AGGG GGAG'):
    adj_list = defaultdict(list)
    for kmer in kmers.split():
        prefix = kmer[:-1]
        sufix = kmer[1:]
        adj_list[prefix].append(sufix)
    return adj_list

def deformat_adj_list(adj_list):
    return {k: v.split() for k, v in [line.split(': ') for line in adj_list.split('\n')]}

def EulerianCycle(adj_list="0: 3\n1: 0\n2: 1 6\n3: 2\n4: 2\n5: 4\n6: 5 8\n7: 9\n8: 7\n9: 6"):
    adj_list = deformat_adj_list(adj_list)
    stack = [k for sublist in adj_list.values() for k in sublist]
    #random.shuffle(stack)
    start = stack.pop()
    cycle = [start]
    while stack:
        while adj_list[cycle[-1]]:
            cycle.append(adj_list[cycle[-1]].pop())
            if cycle[-1] in stack:
                stack.remove(cycle[-1])
        adj_list[cycle[-3]].append(cycle[-2])
        cycle = [cycle[-2]] + cycle[:-2]
    cycle += cycle[:1]
    return cycle

def balance_adj_list(adj_list):
    adj_list = deformat_adj_list(adj_list)
    edges = set([k for sublist in adj_list.values() for k in sublist] + list(adj_list.keys()))
    start, end = None, None
    for e in edges:
        outdegree = len(adj_list[e]) if e in adj_list else 0
        indegree = len([v for v in adj_list.values() if e in v])
        if outdegree < indegree:
            start = e
        if outdegree > indegree:
            end = e
    adj_list[start] = [end] if start not in adj_list else adj_list[start] + [end]
    return adj_list, start, end

def EulerianPath(adj_list="0: 2\n1: 3\n2: 1\n3: 0 4\n6: 3 7\n7: 8\n8: 9\n9: 6"):
    adj_list, s, e = balance_adj_list(adj_list)
    eulerian_cycle = EulerianCycle(format_adj_list(adj_list))[:-1]
    if s and e:
        p1, p2 = ' '.join(eulerian_cycle).split(f'{s} {e}')
        eulerian_path = (e + ' ' + p2 + ' ' + p1 + ' ' + s).split()
    else:
        eulerian_path = ' '.join(eulerian_cycle)
    return eulerian_path

def StringReconstruction(k=4, text='CTTA ACCA TACC GGCT GCTT TTAC'):
    adj_list = DeBruijnGraphFromKmers(text)
    eulerian_path = EulerianPath(format_adj_list(adj_list))
    genome = ''.join([k[0] for k in eulerian_path[:-1]]) + eulerian_path[-1]
    return genome

def kUniversalCircularString(k=4):
    kmers = [''.join(kmer) for kmer in itertools.product('01', repeat=k)]
    #random.shuffle(kmers)
    adj_list = DeBruijnGraphFromKmers(' '.join(kmers))
    eulerian_cycle = EulerianCycle(format_adj_list(adj_list))
    eulerian_cycle = EulerianCycle(format_adj_list(adj_list))[:-(k-1)]
    genome = PathToGenome(' '.join(eulerian_cycle))
    return genome    

def StringSpelledByGappedPatterns(k=4, d=2, text='GACC|GCGC ACCG|CGCC CCGA|GCCG CGAG|CCGG GAGC|CGGA'):
    first_patterns = [pair.split('|')[0] for pair in text.split()]
    second_patterns = [pair.split('|')[1] for pair in text.split()]
    prefix_string = ''.join([p[0] for p in first_patterns]) + first_patterns[-1][1:]
    suffix_string = ''.join([p[0] for p in second_patterns]) + second_patterns[-1][1:]
    for i in range(k + d + 1, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return 1
    return prefix_string + suffix_string[-(k+d):]

def ReconstructFromReadPairs(k=4, d=2, text='GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT'):
    patterns = [tuple(p.split('|')) for p in text.split()]
    adj_list = defaultdict(list)
    for pair in patterns:
        prefix = pair[0][:k-1] + '|' + pair[1][:k-1]
        suffix = pair[0][1:] + '|' + pair[1][1:]
        adj_list[prefix].append(suffix)
    eulerian_path = EulerianPath(format_adj_list(adj_list))
    genome = StringSpelledByGappedPatterns(k, d, ' '.join(eulerian_path))
    return genome

def MaximalNonBranchingPaths(adj_list="1: 2\n2: 3\n3: 4 5\n6: 7\n7: 6"):
    # Paths ← empty list
    # for each node v in Graph
    #     if v is not a 1-in-1-out node
    #         if out(v) > 0
    #             for each outgoing edge (v, w) from v
    #                 NonBranchingPath ← the path consisting of single edge (v, w)
    #                 while w is a 1-in-1-out node
    #                     extend NonBranchingPath by the edge (w, u) 
    #                     w ← u
    #                 add NonBranchingPath to the set Paths
    # for each isolated cycle Cycle in Graph
    #     add Cycle to Paths
    # return Paths
    pass