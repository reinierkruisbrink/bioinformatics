from replication import *
from motifs import *

TASK = 'PatternCount'

if __name__ == "__main__":
    if TASK == 'PatternCount':
        sample_result = PatternCount()
        print(sample_result)
        assert 2 == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = PatternCount(data[0].strip(), data[1].strip())
            print(result)
    elif TASK == 'FrequentWords':
        sample_result = FrequentWords()
        assert 'CATG GCAT' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = FrequentWords(data[0].strip(), int(data[1].strip()))
            print(result)
    elif TASK == 'ReverseComplement':
        sample_result = ReverseComplement()
        assert 'ACCGGGTTTT' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ReverseComplement(data[0].strip())
            print(result)
    elif TASK == 'PatternMatching':
        sample_result = PatternMatching()
        assert '1 3 9' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = PatternMatching(data[0].strip(), data[1].strip())
            #result = PatternMatching('CTTGATCAT', data[0].strip())
            print(result)
    elif TASK == 'ClumpFinding':
        sample_result = ClumpFinding()
        assert 'CGACA GAAGA' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ClumpFinding(data[0].strip(), int(data[1].strip().split()[0]), int(data[1].strip().split()[1]), int(data[1].strip().split()[2]))
            #result = ClumpFinding(data[0].strip(), 9, 500, 3)
            print(result)  
    elif TASK == 'MinSkew':
        sample_result = MinSkew()
        assert '11 24' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = MinSkew(''.join(data).strip())
            print(result)
    elif TASK == 'HammingDistance':
        sample_result = HammingDistance()
        assert 3 == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = str(HammingDistance(data[0].strip(), data[1].strip()))
            print(result)
    elif TASK == 'ApproximatePatternMatching':
        sample_result = ApproximatePatternMatching()
        assert '6 7 26 27' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ApproximatePatternMatching(data[0].strip(), data[1].strip(), int(data[2].strip()))
            print(result)
    elif TASK == 'ApproximatePatternCount':
        sample_result = ApproximatePatternCount()
        assert 4 == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = str(ApproximatePatternCount(data[0].strip(), data[1].strip(), int(data[2].strip())))
            print(result)
    elif TASK == 'Neighbors':
        sample_result = ' '.join(sorted(Neighbors()))
        assert 'AAG ACA ACC ACG ACT AGG ATG CCG GCG TCG' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(Neighbors(data[0].strip(), int(data[1].strip())))
            print(result)
    elif TASK == 'FrequentWordsWithMismatches':
        sample_result = ' '.join(FrequentWordsWithMismatches())
        assert 'ATGT GATG ATGC' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(FrequentWordsWithMismatches(data[0].strip(), int(data[1].strip()[0]), int(data[1].strip()[2])))
            print(result)
    elif TASK == 'FrequentWordsWithMismatchesWithRC':
        sample_result = ' '.join(FrequentWordsWithMismatchesWithRC())
        assert 'ATGT ACAT' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(FrequentWordsWithMismatchesWithRC(data[0].strip(), int(data[1].strip()[0]), int(data[1].strip()[2])))
            print(result)
    elif TASK == 'MotifsEnumeration':
        sample_result = ' '.join(sorted(list(MotifsEnumeration())))
        print(sample_result)
        assert 'ATA ATT GTT TTT' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(MotifsEnumeration(data[1].strip(), int(data[0].split()[0]), int(data[0].split()[1])))
            print(result)
    elif TASK == 'DistanceBetweenPatternAndStrings':
        sample_result = DistanceBetweenPatternAndStrings()
        assert 5 == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = str(DistanceBetweenPatternAndStrings(data[0].strip(), data[1].strip()))
            print(result)
    elif TASK == 'MedianString':
        sample_result = MedianString()
        assert 'GAC' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = MedianString(int(data[0].strip()), data[1].strip())
            print(result)
    elif TASK == 'ProfileMostProbablePattern':
        sample_result = ProfileMostProbablePattern()
        assert 'CCGAG' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ProfileMostProbablePattern(data[0].strip(), int(data[1].strip()), {'A':[float(x) for x in data[2].strip().split()], 'C':[float(x) for x in data[3].strip().split()], 'G':[float(x) for x in data[4].strip().split()], 'T':[float(x) for x in data[5].strip().split()]})
            print(result)        
    elif TASK == 'GreedyMotifSearch':
        sample_result = ' '.join(GreedyMotifSearch())
        assert 'CAG CAG CAA CAA CAA' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(GreedyMotifSearch(data[1].strip(), int(data[0].split()[0]), int(data[0].split()[1])))
            print(result)
    elif TASK == 'GreedyMotifSearchWithPseudocounts':
        sample_result = ' '.join(GreedyMotifSearchWithPseudocounts())
        assert 'TTC ATC TTC ATC TTC' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(GreedyMotifSearchWithPseudocounts(data[1].strip(), int(data[0].split()[0]), int(data[0].split()[1])))
            print(result)
    elif TASK == 'RandomizedMotifSearch':
        sample_result = ' '.join(IterateRandomizedMotifSearch(1000))
        assert 'TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(IterateRandomizedMotifSearch(2000, data[1].strip(), int(data[0].split()[0]), int(data[0].split()[1])))
            print(result)
    elif TASK == 'GibbsSampler':
        #sample_result = ' '.join(IterateGibbsSampler(50))
        #assert 'TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG' == sample_result
        with open(f'data/{TASK}.txt') as f:
            data = f.readlines()
            result = ' '.join(IterateGibbsSampler(50, data[1].strip(), int(data[0].split()[0]), int(data[0].split()[1]), int(data[0].split()[2])))
            print(result)


    else:
        print("Invalid task")
    
    with open('data/output_file.txt', 'w') as f:
        f.write(result)

    # print(HammingDistance('TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC', 'GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA'))
    # print(MinSkew('CATTCCAGTACTTCGATGATGGCGTGAAGA'))
    # print(ApproximatePatternCount('TGT', 'CGTGACAGTGTATGGGCATCTTT', 1))
    # print(len(Neighbors('ACGT', 2)))
    # print(MedianString(7, 'CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCC CGGCGCTAATCCTAGTGCCC GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACC ACCACGGGTGGCTAGTTTC GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAAT CGGCTTCAAATCCTACACAG'))

    # dna = 'ATGAGGTC GCCCTAGA AAATAGAT TGGTGCTA'.split()
    # bm = 'GTC CCC ATA GCT'.split()
    # print(' '.join(Motifs(Profile(bm), dna, 3)))