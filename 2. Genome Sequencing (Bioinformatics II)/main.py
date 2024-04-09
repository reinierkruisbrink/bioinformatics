from dna_sequencing import *
from peptides_sequencing import *

TASK = 'CyclopeptideSequencing'

if __name__ == "__main__":
    if TASK == 'Composition':
        sample_result = ' '.join(composition())
        assert 'CAATC AATCC ATCCA TCCAA CCAAC' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(composition(int(data[0]), data[1].strip()))
            print(result)
    elif TASK == 'PathToGenome':
        sample_result = ''.join(PathToGenome())
        assert 'ACCGAAGCT' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ''.join(PathToGenome(data[0]))
            print(result)
    elif TASK == 'OverlapGraph':
        sample_result = format_adj_list(OverlapGraph())
        assert """GCATG: CATGC\nCATGC: ATGCG\nAGGCA: GGCAT GGCAC\nGGCAT: GCATG""" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = format_adj_list(OverlapGraph(data[0].strip()))
            print(result)
    elif TASK == 'DeBruijnGraph':
        sample_result = format_adj_list(DeBruijnGraph())
        print(sample_result)
        assert """AAG: AGA AGA\nAGA: GAT\nGAT: ATT\nATT: TTC\nTTC: TCT\nTCT: CTC CTA\nCTC: TCT\nCTA: TAA\nTAA: AAG""" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = format_adj_list(DeBruijnGraph(int(data[0]), data[1].strip()))
            print(result)
    elif TASK == 'DeBruijnGraphFromKmers':
        sample_result = format_adj_list(DeBruijnGraphFromKmers())
        print(sample_result)
        assert """GAG: AGG\nCAG: AGG AGG\nGGG: GGG GGA\nAGG: GGG\nGGA: GAG""" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = format_adj_list(DeBruijnGraphFromKmers(data[0].strip()))
            print(result)
    elif TASK == 'EulerianCycle':
        sample_result = ' '.join(EulerianCycle())
        #assert "6 8 7 9 6 5 4 2 1 0 3 2 6" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(EulerianCycle(''.join(data).strip()))
            print(result)
    elif TASK == 'EulerianPath':
        sample_result = ' '.join(EulerianPath())
        #assert "6 7 8 9 6 3 0 2 1 3 4" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(EulerianPath(''.join(data).strip()))
            print(result)
    elif TASK == 'StringReconstruction':
        sample_result = StringReconstruction()
        assert "GGCTTACCA" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = StringReconstruction(int(data[0].strip()), ''.join(data[1:]))
            print(result)
    elif TASK == 'kUniversalCircularString':
        sample_result = kUniversalCircularString()
        #assert "00111010" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = kUniversalCircularString(int(data[0].strip()))
            print(result)
    elif TASK == 'StringSpelledByGappedPatterns':
        sample_result = StringSpelledByGappedPatterns()
        assert "GACCGAGCGCCGGA" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = StringSpelledByGappedPatterns(int(data[0].split()[0]), int(data[0].split()[1]), data[1].strip())
            print(result)
    elif TASK == 'ReconstructFromReadPairs':
        sample_result = ReconstructFromReadPairs()
        assert "GTGGTCGTGAGATGTTGA" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ReconstructFromReadPairs(int(data[0].split()[0]), int(data[0].split()[1]), data[1].strip())
            print(result)
    elif TASK == 'MaximalNonBranchingPaths':
        sample_result = MaximalNonBranchingPaths()
        assert "1 2 3\n3 4\n3 5\n7 6 7" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = MaximalNonBranchingPaths(''.join(data).strip())
            print(result)
    elif TASK == 'ProteinTranslation':
        sample_result = ProteinTranslation()[:-1]
        assert "MAMAPRTEINSTRING" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ProteinTranslation(data[0].strip())[:-1]
            print(result)
    elif TASK == 'PeptideEncoding':
        sample_result = '\n'.join(PeptideEncoding())
        assert "ATGGCC\nGGCCAT\nATGGCC" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = '\n'.join(PeptideEncoding(data[0].strip(), data[1].strip()))
            #result = len(PeptideEncoding(data[0].strip(), 'VKLFPWFNQY'))
            print(result)
    elif TASK == 'CyclopeptideSequencingProblem':
        sample_result = CyclopeptideSequencingProblem()
        assert "980597910" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = CyclopeptideSequencingProblem(int(data[0].strip()))
            print(result)
    elif TASK == 'LinearSpectrum':
        sample_result = ' '.join([str(i) for i in LinearSpectrum()])
        assert "0 113 114 128 129 242 242 257 370 371 484" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join([str(i) for i in LinearSpectrum(data[0].strip())])
            print(result)
    elif TASK == 'CyclicSpectrum':
        sample_result = ' '.join([str(i) for i in CyclicSpectrum()])
        assert "0 113 114 128 129 227 242 242 257 355 356 370 371 484" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join([str(i) for i in CyclicSpectrum(data[0].strip())])
            print(result)
    elif TASK == 'CountingPeptides':
        sample_result = str(CountingPeptides())
        assert "14712706211" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = str(CountingPeptides(int(data[0].strip())))
            print(result)
    elif TASK == 'CountingSubpeptides':
        sample_result = str(CountingSubpeptides())
        assert "11" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = str(CountingSubpeptides(int(data[0].strip())))
            print(result)
    elif TASK == 'CyclopeptideSequencing':
        sample_result = ' '.join(CyclopeptideSequencing())
        assert "186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(CyclopeptideSequencing(data[0].strip()))
            print(result)

    else:
        print("Invalid task")

    with open('data/output_file.txt', 'w') as f:
        f.write(result)

