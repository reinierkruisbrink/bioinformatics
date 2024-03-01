from collections import Counter, defaultdict


def PatternCount(sequence='GCGCG', Pattern='GCG'):
    count = 0
    for i in range(len(sequence)-len(Pattern)+1):
        if sequence[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def FrequentWords(sequence='ACGTTGCATGTCGCATGATGCATGAGAGCT', k=4):
    words = []
    freq = dict(Counter([sequence[i:i+k] for i in range(len(sequence)) if len(sequence[i:i+k])==k]))
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return ' '.join(words[::-1])

def ReverseComplement(sequence='AAAACCCGGT'):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([complement[i] for i in sequence[::-1]])

def PatternMatching(Pattern='ATAT', sequence='GATATATGCATATACTT'):
    return ' '.join([str(i) for i in range(len(sequence)-len(Pattern)+1) if sequence[i:i+len(Pattern)] == Pattern])

def ClumpFinding(sequence='CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', k=5, L=50, t=4):
    return ' '.join([key for key, value in dict(Counter([sequence[i:i+k] for i in range(len(sequence)) if len(sequence[i:i+k])==k])).items() if value >= t])

def MinSkew(sequence='TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'):
    skew = [0]
    for i in range(len(sequence)):
        if sequence[i] == 'C':
            skew.append(skew[i]-1)
        elif sequence[i] == 'G':
            skew.append(skew[i]+1)
        else:
            skew.append(skew[i])
    m = min(skew)
    return ' '.join([str(i) for i in range(len(skew)) if skew[i] == m])

def HammingDistance(sequence1='GGGCCGTTGGT', sequence2='GGACCGTTGAC'):
    dist = sum([1 for a,b in zip(sequence1, sequence2) if a != b])
    return dist

def ApproximatePatternMatching(Pattern='ATTCTGGA', sequence='CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', d=3):
    positions = []
    for i in range(len(sequence)-len(Pattern)+1):
        if HammingDistance(Pattern, sequence[i:i+len(Pattern)]) <= d:
            positions.append(i)
    return ' '.join([str(p) for p in positions])

def ApproximatePatternCount(Pattern='GAGG', sequence='TTTAGAGCCTTCAGAGG', d=2):
    count = len(ApproximatePatternMatching(Pattern, sequence, d).split())
    return count

def Neighbors(Pattern='ACG', d=1):
    if d == 0:
        return Pattern
    if len(Pattern) == 1:
        return ['A', 'C', 'G', 'T']
    neighborhood = []
    suffixNeighbors = Neighbors(Pattern[1:], d)
    for seq in suffixNeighbors:
        if HammingDistance(Pattern[1:], seq) < d:
            for n in ['A', 'C', 'G', 'T']:
                neighborhood.append(n+seq)
        else:
            neighborhood.append(Pattern[0]+seq)
    return neighborhood

def FrequentWordsWithMismatches(sequence='ACGTTGCATGTCGCATGATGCATGAGAGCT', k=4, d=1):
    Patterns = []
    freqMap = defaultdict(int)
    for i in range(len(sequence)-k+1):
        neighborhood = Neighbors(sequence[i:i+k], d)
        for neighbor in neighborhood:
            freqMap[neighbor] += 1
    m = max(freqMap.values())
    for Pattern in freqMap:
        if freqMap[Pattern] == m:
            Patterns.append(Pattern)
    return Patterns

def FrequentWordsWithMismatchesWithRC(sequence='ACGTTGCATGTCGCATGATGCATGAGAGCT', k=4, d=1):
    Patterns = []
    freqMap = {}
    for i in range(len(sequence)-k+1):
        neighborhood = Neighbors(sequence[i:i+k], d)
        for neighbor in neighborhood:
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1

    for Pattern in freqMap:
        if ReverseComplement(Pattern) in freqMap:
            freqMap[Pattern] += freqMap[ReverseComplement(Pattern)]

    m = max(freqMap.values())
    for Pattern in freqMap:
        if freqMap[Pattern] == m:
            Patterns.append(Pattern)
            Patterns.append(ReverseComplement(Pattern))
    return Patterns
