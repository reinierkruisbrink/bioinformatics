from sequencing import *

TASK = 'DeBruijnGraphFromKmers'

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
        sample_result = OverlapGraph()
        assert """GCATG: CATGC
CATGC: ATGCG
AGGCA: GGCAT GGCAC
GGCAT: GCATG""" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = OverlapGraph(data[0].strip())
            print(result)
    elif TASK == 'DeBruijnGraph':
        sample_result = DeBruijnGraph()
        print(sample_result)
        assert """AAG: AGA AGA
AGA: GAT
GAT: ATT
ATT: TTC
TTC: TCT
TCT: CTC CTA
CTC: TCT
CTA: TAA
TAA: AAG""" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = DeBruijnGraph(int(data[0]), data[1].strip())
            print(result)
    elif TASK == 'DeBruijnGraphFromKmers':
        sample_result = DeBruijnGraphFromKmers()
        print(sample_result)
        assert """GAG: AGG
CAG: AGG AGG
GGG: GGG GGA
AGG: GGG
GGA: GAG""" == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = DeBruijnGraphFromKmers(data[0].strip())
            print(result)

    else:
        print("Invalid task")
    
    with open('data/output_file.txt', 'w') as f:
        f.write(result)