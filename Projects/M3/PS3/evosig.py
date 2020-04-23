## PS3 problem 2: Discovery of evolutionary signatures.

## Given all intergeneic regions of a piece of DNA, and which nucleotides
## are conserved, compute the frequency and and conservation of k-mers
## (k = 6 for part A). 

## Output for part A:
##  -A list of the 50 most frequent k-mers
##  -A list of the 50 most conserved k-mers (we will say all nucleotides in
##      the k-mer must be conserved for the k-mer to count as conserved)


import operator

codons = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'n', 'TAG':'n', 
        'TGC':'C', 'TGT':'C', 'TGA':'n', 'TGG':'W', 
    } 


def readdata(file):

    data = []
    for line in open(file,'r'):
        if line[-1] == "#":
            line = line[0:-1]
        data.append(line)
    return data


def getKmers(kmerlen, seq, conserved):

    kmers = {}

    for i in range(len(seq) - kmerlen + 1):

        kmer = seq[i:i+kmerlen]

        # if kmer contains non-base characters, replace with 'n'
        # (will be useful for part B)
        # if "-" in kmer or "#" in kmer:
        #     kmer.replace("-", "n")
        #     kmer.replace("#", "n")

        conservedSeq = conserved[i:i+kmerlen]
        numConserved = conservedSeq.count("*")

        if kmer in kmers:
            kmers[kmer][0].append(i)
            kmers[kmer][1].append(numConserved)
        else:
            kmers[kmer] = [[i], [numConserved]]
    
    return kmers


def nucleotideToAminoAcid(seq):

    N = len(seq) / 3
    aaByPosition = {}
    bases = ['A','C','G','T']

    for i in range(N):

        cod = seq[3*i : 3*i+3]

        if "-" not in cod and "#" not in cod:

            aa = codons[cod]
            aaByPosition.setdefault(i, []).append(aa)

        else:

            # find possible codons
            possibleCodons = []
            for j in range(len(bases)):
                if "#" in cod:
                    possibleCodons.append(cod.replace("#", bases[j]))
                if "-" in cod:
                    possibleCodons.append(cod.replace("-", bases[j]))

            # get possible aa
            possibleAA = []
            for c in possibleCodons:
                possibleAA.append(codons[c])

            # remove duplicates and add to dictionary
            aaByPosition[i] = list(dict.fromkeys(possibleAA))

    # enumerate all possible sequences

    numAASeqs = 1
    for p in aaByPosition:
        numAASeqs *= len(aaByPosition[p])

    
    aaSeqs = ["" for i in range(numAASeqs)]
    interval = numAASeqs
    for i in range(len(aaByPosition)):

        numPossibilities = len(aaByPosition[i])
        interval /= numPossibilities

        for j in range(numAASeqs):
            aaSeqs[j] = aaSeqs[j] + aaByPosition[i][(j/interval) % numPossibilities]
    
    # remove duplicates and add to dictionary
    aaSeqs = list(dict.fromkeys(aaSeqs))

    return aaSeqs


# for part B, find most conserved k-mers in yeast_motifs.txt
def findKmersInYeast(kmers):

    # read yeast data from file
    yeastfile = "yeast_motifs.txt"
    yeastdata = {}
    for line in open(yeastfile,'r'):
        if line[-1] == "\n":
            line = line[0:-1]
        motifdata = line.split(" ")
        if len(motifdata) < 2:
            continue
        yeastdata[motifdata[1]] = motifdata[0]
    
    print "Kmers found in yeast:"
    commonKmers = {}
    for kmer in kmers:
        # convert kmers from bases to amino acids
        aaSeqs = nucleotideToAminoAcid(kmer)
        for kmer in aaSeqs:
            # check for kmers in yeast
            if kmer in yeastdata:
                commonKmers[kmer] = yeastdata[kmer]
                print kmer, " : ", yeastdata[kmer]
    
    return commonKmers
    

def main():

    ## READ DATA FROM FILE

    intergeneicfile = "allinter"
    annotationfile = "allintercons"

    intergenicList = readdata(intergeneicfile)
    annotationList = readdata(annotationfile)

    if len(intergenicList) != len(annotationList):
        print "Number of annotations is not equal to the number of sequences"
        exit(1)

    # since we know it's only one sequence, just take first element in list
    intergenic = intergenicList[0]
    annotation = annotationList[0]
    kmerlen = 6

    ## EXTRACT KMERS FROM DATA

    kmers = getKmers(kmerlen, intergenic, annotation)
    # kmers = { kmer : [ [indexLocations] , [conservationCounts] ] }

    ## FIND THE TOP 50 FREQUENT AND TOP 50 CONSERVED KMERS

    mostFrequent = {}
    mostConserved = {}

    N = 50  # top entries

    for kmer in kmers:

        kmerCount = len(kmers[kmer][0])
        kmerProporitonConserved = \
            float(sum(kmers[kmer][1])) / (kmerlen*len(kmers[kmer][1]))  # proportion of conserved instances

        # if list is less than 50, add
        if len(mostFrequent) < N:
            mostFrequent[kmer] = kmerCount

        # else, remove one of the min entries and add
        else:
            minCount = min(mostFrequent.values())
            if kmerCount > minCount:
                minKmers = [key for key in mostFrequent if mostFrequent[key] == minCount]
                mostFrequent.pop(minKmers[0])
                mostFrequent[kmer] = kmerCount
        
        # do the same for conserved
        if len(mostConserved) < N:
            mostConserved[kmer] = kmerProporitonConserved
        else:
            minCount = min(mostConserved.values())
            if kmerProporitonConserved > minCount:
                minKmers = [key for key in mostConserved if mostConserved[key] == minCount]
                mostConserved.pop(minKmers[0])
                mostConserved[kmer] = kmerProporitonConserved
    
    ## SORT TOP 50 DICTIONARIES

    mostFrequentList = sorted(mostFrequent.items(), key=operator.itemgetter(1), reverse=True)
    mostConservedList = sorted(mostConserved.items(), key=operator.itemgetter(1), reverse=True)

    ## PRINT RESULTS

    print "Top 50 most FREQUENT:"
    for (kmer,count) in mostFrequentList:
        print kmer, " : ", count 
    print ""
    print "Top 50 most CONSERVED:"
    for (kmer,count) in mostConservedList:
        print kmer, " : ", count

    ## PART B, CHECK FOR MOST CONSERVED K-MERS IN YEAST

    # NOTE: it's not clear if these yeast file contains amino 
    # sequences or nucleotide sequences. I assumed it was nucleotides
    # since there are letters other than the 4 bases. Wrote a function
    # to convert codons to AA.

    yeastKmers = findKmersInYeast(mostConserved)


if __name__ == "__main__":
    main()
