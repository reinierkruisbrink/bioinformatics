import itertools
from collections import Counter
from Bio.Seq import Seq


# Codon to Amino Acid
C2A = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T',  'ACC': 'T', 'ACG': 'T', 'ACU': 'T', 'AGA': 'R', 'AGC': 'S',  'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M',  'AUU': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H',  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P', 'CGA': 'R',  'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L', 'CUC': 'L',  'CUG': 'L', 'CUU': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',  'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V',  'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'UAA': '*', 'UAC': 'Y',  'UAG': '*', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S',  'UCU': 'S', 'UGA': '*', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C',  'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}
# Amino Acid to Codon
A2C = {'K': ['AAA', 'AAG'], 'N': ['AAC', 'AAU'], 'T': ['ACA', 'ACC', 'ACG', 'ACU'], 'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'], 'S': ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'], 'I': ['AUA', 'AUC', 'AUU'], 'M': ['AUG'], 'Q': ['CAA', 'CAG'], 'H': ['CAC', 'CAU'], 'P': ['CCA', 'CCC', 'CCG', 'CCU'], 'L': ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'], 'E': ['GAA', 'GAG'], 'D': ['GAC', 'GAU'], 'A': ['GCA', 'GCC', 'GCG', 'GCU'], 'G': ['GGA', 'GGC', 'GGG', 'GGU'], 'V': ['GUA', 'GUC', 'GUG', 'GUU'], '*': ['UAA', 'UAG', 'UGA'], 'Y': ['UAC', 'UAU'], 'C': ['UGC', 'UGU'], 'W': ['UGG'], 'F': ['UUC', 'UUU']}

# Integer Mass Table
A2M = {'A': 71, 'C': 103, 'D': 115, 'E': 129, 'F': 147, 'G': 57, 'H': 137, 'I': 113, 'K': 128, 'L': 113, 'M': 131, 'N': 114, 'P': 97, 'Q': 128, 'R': 156, 'S': 87, 'T': 101, 'V': 99, 'W': 186, 'Y': 163}
# Mass to Amino Acid
M2A = {71: 'A', 103: 'C', 115: 'D', 129: 'E', 147: 'F', 57: 'G', 137: 'H', 113: 'I', 128: 'K', 131: 'M', 114: 'N', 97: 'P', 156: 'R', 87: 'S', 101: 'T', 99: 'V', 186: 'W', 163: 'Y'}

def ProteinTranslation(rna='AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'):
    peptide = ''
    for i in range(0, len(rna), 3):
        peptide += C2A[rna[i:i+3]]
        if peptide[-1] == '*':
            break
    return peptide

def PeptideEncoding(dna='ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', peptide='MA'):
    results = []
    rna = dna.replace('T', 'U')
    for i in range(len(dna) - len(peptide)*3+1):
        codon = Seq(rna[i:i+len(peptide)*3])
        if ProteinTranslation(codon) == peptide:
            results.append(dna[i:i+len(peptide)*3])
        if ProteinTranslation(codon.reverse_complement_rna()) == peptide:
            results.append(dna[i:i+len(peptide)*3])
    return results

def CyclopeptideSequencingProblem(n=31315):
    return str(n**2-n)

def LinearSpectrum(peptide='NQEL', mass_table=A2M):
    prefix_mass = [0]
    for i in range(1, len(peptide)+1):
        for s in mass_table:
            if s == peptide[i-1]:
                prefix_mass.append(prefix_mass[i-1] + mass_table[s])
    linear_spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i+1, len(peptide)+1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)

def CyclicSpectrum(peptide='LEQN', mass_table=A2M):
    prefix_mass = [0]
    for i in range(1, len(peptide)+1):
        for s in mass_table:
            if s == peptide[i-1]:
                prefix_mass.append(prefix_mass[i-1] + mass_table[s])
    peptide_mass = prefix_mass[-1]
    cyclic_spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i+1, len(peptide)+1):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cyclic_spectrum)

def CountingPeptides(n=1024):
    mass_list = list(set(A2M.values()))
    table = [1] + [0]*n

    def count_peptides(m):
        if m < 0:
            return 0
        if table[m] != 0:
            return table[m]
        tmp = 0
        for mass in mass_list:
            tmp += count_peptides(m - mass)
        table[m] = tmp
        return tmp

    return count_peptides(n)

def CountingSubpeptides(n=4):
    return sum([i for i in range(n+1)]) + 1

def Expand(peptides, extended=False):
    if extended:
        mass_list = list(range(57, 201))
    else:
        mass_list = list(set(A2M.values()))
    expanded_peptides = set()
    for p in peptides:
        if p == '':
            for m in mass_list:
                expanded_peptides.add(str(m))
        else:
            for m in mass_list:
                expanded_peptides.add(p+ '-' + str(m))
    return expanded_peptides

def CycloSpectrum(aa_list):
    # Same as CyclicSpectrum but for mass list instead of peptide string
    prefix_mass = [0]
    for i in range(len(aa_list)):
        prefix_mass.append(prefix_mass[i] + aa_list[i])
    cyclic_spectrum = [0]
    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)+1):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(aa_list):
                cyclic_spectrum.append(prefix_mass[-1]-(prefix_mass[j]-prefix_mass[i]))
    return Counter(cyclic_spectrum)

def LineoSpectrum(aa_list):
    # Same as LinearSpectrum but for mass list instead of peptide string
    prefix_mass = [0]
    for i in range(len(aa_list)):
        prefix_mass.append(prefix_mass[i] + aa_list[i])
    linear_spectrum = [0]
    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)+1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return Counter(linear_spectrum)

def IsConsistent(aa_list, spectrum_dict):
    candidate_spectrum_dict = LineoSpectrum(aa_list)
    for key, value in candidate_spectrum_dict.items():
        if value > spectrum_dict.get(key, 0):
            return False
    return True

def CyclopeptideSequencing(spectrum='0 113 128 186 241 299 314 427'):
    spectrum = [int(i) for i in spectrum.split()]
    spectrum_dict = Counter(spectrum)
    parent_mass = max(spectrum)
    candidate_peptides = {''}
    final_peptides = []
    while len(candidate_peptides) > 0:
        candidate_peptides = Expand(candidate_peptides)
        deletions = []
        for peptide in candidate_peptides:
            aa_list = [int(aa) for aa in peptide.split('-')]
            if sum(aa_list) == parent_mass:
                if CycloSpectrum(aa_list) == spectrum_dict and peptide not in final_peptides:
                    final_peptides.append(peptide)
                deletions.append(peptide)
            elif not IsConsistent(aa_list, spectrum_dict):
                deletions.append(peptide)
        for p in deletions:
                candidate_peptides.remove(p)
    return sorted(final_peptides)[::-1]

def CyclopeptideScoring(peptide='NQEL', spectrum='0 99 113 114 128 227 257 299 355 356 370 371 484'):
    if peptide == '':
        return 0
    spectrum = [int(i) for i in spectrum.split()]
    spectrum_dict = Counter(spectrum)
    # aa_list = [A2M[aa] for aa in peptide] #peptide string
    aa_list = [int(aa) for aa in peptide.split('-')] #mass list
    candidate_spectrum_dict = CycloSpectrum(aa_list)
    score = 0
    for key, value in candidate_spectrum_dict.items():
        score += min(value, spectrum_dict.get(key, 0))
    return score

def LinearScore(peptide='NQEL', spectrum='0 99 113 114 128 227 257 299 355 356 370 371 484'):
    spectrum = [int(i) for i in spectrum.split()]
    spectrum_dict = Counter(spectrum)
    # aa_list = [A2M[aa] for aa in peptide] #peptide string
    aa_list = [int(aa) for aa in peptide.split('-')] #mass list
    candidate_spectrum_dict = LineoSpectrum(aa_list)
    score = 0
    for key, value in candidate_spectrum_dict.items():
        score += min(value, spectrum_dict.get(key, 0))
    return score

def Trim(leaderboard='LAST ALST TLLT TQAS', spectrum='0 71 87 101 113 158 184 188 259 271 372', n=2):
    linear_scores = [(peptide, LinearScore(peptide, spectrum)) for peptide in leaderboard.split()]
    linear_scores.sort(key=lambda x: x[1], reverse=True)
    return set([peptide for peptide, _ in linear_scores[:n]])

def LeaderboardCyclopeptideSequencing(n=10, spectrum='0 71 113 129 147 200 218 260 313 331 347 389 460', extended=False, parent_mass=None):
    if parent_mass is None:
        parent_mass = max([int(i) for i in spectrum.split()])
    leaderboard = {''}
    leader_peptide = ['']
    best_score = 0
    while len(leaderboard) > 0:
        leaderboard = Expand(leaderboard, extended=extended)
        deletions = []
        for peptide in leaderboard:
            aa_list = [int(aa) for aa in peptide.split('-')]
            if sum(aa_list) == parent_mass:
                peptide_score = CyclopeptideScoring(peptide, spectrum)
                if peptide_score > best_score:
                    leader_peptide = [peptide]
                    best_score = peptide_score
                elif peptide_score == best_score:
                    leader_peptide.append(peptide)
            elif sum(aa_list) > parent_mass:
                deletions.append(peptide)
        for peptide in deletions: 
            leaderboard.remove(peptide)
        leaderboard = Trim(' '.join(leaderboard), spectrum, n)
    return leader_peptide

def SpectralConvolution(spectrum='0 137 186 323'):
    spectrum = [int(i) for i in spectrum.split()]
    convolution = [str(abs(m1-m2)) for m1, m2 in itertools.combinations(spectrum, 2) if m1-m2 != 0]
    return convolution

def ConvolutionCyclopeptideSequencing(m=20, n=60, spectrum='57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493'):
    parent_mass = max([int(i) for i in spectrum.split()])
    spectral_convolution = sorted(Counter(SpectralConvolution(spectrum)).items(), key=lambda x: x[1], reverse=True)
    print(spectral_convolution)
    min_multiplicity = spectral_convolution[m-1][1]
    spectrum = ''
    for mass, multiplicity in spectral_convolution:
        if multiplicity >= min_multiplicity and int(mass) >= 57 and int(mass) < 200:
            spectrum += (mass + ' ') * multiplicity
    print(spectrum)
    return LeaderboardCyclopeptideSequencing(n, spectrum, extended=True, parent_mass=parent_mass)
