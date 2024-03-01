import math
import random
import itertools
from replication import *

def MotifsEnumeration(dna='ATTTGGC TGCCTTA CGGTATC GAAAATT', k=3, d=1):
    dna = dna.split()
    patterns = set()
    for i in range(len(dna[0]) - k + 1):
        pattern = dna[0][i:i+k]
        neighborhood = Neighbors(pattern, d)
        for neighbor in neighborhood:
            if all([ApproximatePatternMatching(neighbor, dna[j], d) for j in range(1, len(dna))]):
                patterns.add(neighbor)
    return patterns

def H(p):
    return -sum([p_i * math.log(p_i, 2) for p_i in p])

def DistanceBetweenPatternAndStrings(Pattern='AAA', Dna='TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT'):
    k = len(Pattern)
    distance = 0
    for Text in Dna.split():
        Distance = float('inf')
        for i in range(len(Text) - k + 1):
            PatternPrime = Text[i:i+k]
            if Distance > HammingDistance(Pattern, PatternPrime):
                Distance = HammingDistance(Pattern, PatternPrime)
        distance += Distance
    return distance

def MedianString(k=3, dna='AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTTCGGGACAG'):
    distance = float('inf')
    Median = ''
    for Pattern in [''.join(base) for base in itertools.product('ACGT', repeat=k)]:
        d = DistanceBetweenPatternAndStrings(Pattern, dna)
        if distance > d:
            distance = d
            Median = Pattern
    return Median

def Pr(sequence, Profile):
    p = 1.0
    for i in range(len(sequence)):
        p *= Profile[sequence[i]][i]
    return p

def ProfileMostProbablePattern(sequence='ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', k=5, profile={'A':[.2,.2,.3,.2,.3], 'C':[.4,.3,.1,.5,.1], 'G':[.3,.3,.5,.2,.4], 'T':[.1,.2,.1,.1,.2]}):
    max_probability = -1.0
    most_probable_kmer = ""
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        probability = Pr(kmer, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer

def Profile(Motifs):
    k = len(Motifs[0])
    profile = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            base = Motifs[i][j]
            profile[base][j] += 1/t
    return profile

def Count(Motifs):
    k = len(Motifs[0])
    count = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            base = Motifs[i][j]
            count[base][j] += 1
    return count

def Consensus(Motifs, pseudocounts=False):
    k = len(Motifs[0])
    count = Count(Motifs) if not pseudocounts else CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentBase = ""
        for base in "ACGT":
            if count[base][j] > m:
                m = count[base][j]
                frequentBase = base
        consensus += frequentBase
    return consensus

def Score(Motifs, pseudocounts=False):
    cons = Consensus(Motifs, pseudocounts)
    return sum([HammingDistance(cons, seq) for seq in Motifs])

def GreedyMotifSearch(Dna='GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG', k=3, t=5):
    Dna = Dna.split()
    n = len(Dna[0])
    BestMotifs = [seq[:k] for seq in Dna]
    for i in range(n - k + 1):
        Motifs = [Dna[0][i:i+k]]
        for j in range(1, t):
            profile = Profile(Motifs)
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, profile))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def ProfileWithPseudocounts(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)
    profile = {'A':[1/(t+4)]*k, 'C':[1/(t+4)]*k, 'G':[1/(t+4)]*k, 'T':[1/(t+4)]*k}
    for i in range(t):
        for j in range(k):
            base = Motifs[i][j]
            profile[base][j] += 1/(t+4)
    return profile

def CountWithPseudocounts(Motifs):
    k = len(Motifs[0])
    count = {'A':[1]*k, 'C':[1]*k, 'G':[1]*k, 'T':[1]*k}
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            base = Motifs[i][j]
            count[base][j] += 1
    return count

def GreedyMotifSearchWithPseudocounts(Dna='GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG', k=3, t=5):
    Dna = Dna.split()
    n = len(Dna[0])
    BestMotifs = [seq[:k] for seq in Dna]
    for i in range(n-k+1):
        Motifs = [Dna[0][i:i+k]]
        for j in range(1, t):
            profile = ProfileWithPseudocounts(Motifs)
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, profile))
        if Score(Motifs, pseudocounts=True) < Score(BestMotifs, pseudocounts=True):
            BestMotifs = Motifs
    return BestMotifs

def Motifs(Profile, Dna, k):
    Motifs = []
    for seq in Dna:
        Motifs.append(ProfileMostProbablePattern(seq, k, Profile))
    return Motifs

def RandomMotifs(Dna, k):
    Motifs = []
    for seq in Dna:
        r = random.randint(1, len(Dna[0])-k)
        Motifs.append(seq[r:r+k])
    return Motifs

def RandomizedMotifSearch(Dna='CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA', k=8, t=5):
    Dna = Dna.split()
    M = RandomMotifs(Dna, k)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna, k)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs, Score(BestMotifs)
        
def IterateRandomizedMotifSearch(N, Dna='CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA', k=8, t=5):
    BestScore = float('inf')
    for _ in range(N):
        Motifs, Score = RandomizedMotifSearch(Dna, k, t)
        if Score < BestScore:
            BestMotifs = Motifs
            BestScore = Score
    return BestMotifs

def GibbsSampler(Dna='CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA', k=8, t=5, N=1000):
    Dna = Dna.split()
    BestMotifs = []
    Motifs = RandomMotifs(Dna, k)
    BestMotifs = Motifs
    for _ in range(N):
        i = random.randint(0, t-1)
        Profile = ProfileWithPseudocounts(Motifs[:i] + Motifs[i+1:])
        Motifs[i] = ProfileMostProbablePattern(Dna[i], k, Profile)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs, Score(BestMotifs)

def IterateGibbsSampler(I, Dna='CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA', k=8, t=5, N=1000):
    BestScore = float('inf')
    for _ in range(I):
        Motifs, Score = GibbsSampler(Dna, k, t, N)
        if Score < BestScore:
            BestMotifs = Motifs
            BestScore = Score
    return BestMotifs
