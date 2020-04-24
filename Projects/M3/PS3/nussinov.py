## PS3 problem 3: RNA folding, Nussinov algorithm

# Reference: https://math.mit.edu/classes/18.417/Slides/rna-prediction-nussinov.pdf

import numpy as np
import matplotlib.pyplot as plt

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'U' : 3}

# base pair scoring matrix
S = [
    #   A   G   C   U
    [   0,  0,  0,  1],    # A
    [   0,  0,  1,  0],    # G
    [   0,  1,  0,  0],    # C
    [   1,  0,  0,  0]     # U
]

# store scores in F
# F = [
#     1 2 . . . N
#   [. . . . . . .],    1
#   [. . . . . . .],    2
#   . . .
#   [. . . . . . .]     N
# ]


def traceback(F, seq):

    i = 0
    j = len(seq)-1

    structure = []

    while i <= j:

        # set diagonal (left, down) as default
        iNext = i+1
        jNext = j-1
        maxValue = F[iNext][jNext]

        # if left is higher, go left
        if F[i][j-1] > maxValue:
            iNext = i
            jNext = j-1
            maxValue = F[iNext][jNext]

        # if down is higher, go down
        if F[i+1][j] > maxValue:
            iNext = i+1
            jNext = j
            maxValue = F[iNext][jNext]
        
        # only add a pair if they went diagonally
        if iNext == i+1 and jNext == j-1:
            structure.append((i,j))

        i = iNext
        j = jNext
    
    return structure


def calculateScore(i, j, seq, F):

    Si = seq[i]
    Sj = seq[j]
    score1 = F[i+1][j-1] + S[base_idx[Si]][base_idx[Sj]]
    score2 = [F[i][k] + F[k+1][j] for k in range(i, j)]
    return max(score1, max(score2))


def Nussinov(seq):

    N = len(seq)    
    F = [ [0.0 for i in range(N)] for j in range(N)]

    for l in range(1, N):
        for i in range(0, N-l):
            j = i + l
            F[i][j] = calculateScore(i, j, seq, F)
    
    score = F[0][N-1]
    
    return F, score


def GenerateRandomRna(n, gcContent=0.5):

    rna = ""
    for i in range(n):
        r = np.random.uniform()
        if r < gcContent/2:
            rna += 'C'
        elif r < gcContent:
            rna += 'G'
        elif r < (1.0-gcContent/2):
            rna += 'A'
        else:
            rna += 'U'

    return rna


def printResults(F, score, structure):

    print "Score = ", score
    print "F = "
    for i in range(len(F)):
        print F[i]
    print ""

    print "Set of binding pairs: ", structure
    print "-" * 30


##### FUNCTIONS TO ANSWER THE DIFFERENT PARTS OF PROBLEM 3 #####

# part A
def runSingleSample(n=100):

    ## GENERATE RANDOM RNAS BY UNIFORM DISTRIBUTION

    rna = GenerateRandomRna(n)

    ## RUN NUSSINOV ALGORITHM

    F, score = Nussinov(rnas100[i])

    ## TRACEBACK TO GET SET OF BINDING PAIRS

    S = traceback(F, rna)

    ## PRINT RESULTS (for single sequence only)
    
    printResults(F, score, S)


# part B
def runLen100():

    ## GENERATE RANDOM RNAS BY UNIFORM DISTRIBUTION

    numSamples = 1000
    rnas100 = []
    for i in range(numSamples):
        rnas100.append(GenerateRandomRna(100))

    ## RUN NUSSINOV ALGORITHM

    scoreSum = 0.0
    for i in range(numSamples):
        F, score = Nussinov(rnas100[i])
        scoreSum += score

    averageScore = scoreSum/numSamples

    print "Average score for length 100 = ", averageScore


# part C
def runVariousLengths():

    ## GENERATE RANDOM RNAS BY UNIFORM DISTRIBUTION

    numSamples = 50
    rnas10 = []
    rnas50 = []
    rnas100 = []
    for i in range(numSamples):
        rnas10.append(GenerateRandomRna(10))
        rnas50.append(GenerateRandomRna(50))
        rnas100.append(GenerateRandomRna(100))

    ## RUN NUSSINOV ALGORITHM

    sums = [0.0, 0.0, 0.0]
    for i in range(numSamples):
        F10, score10 = Nussinov(rnas10[i])
        F50, score50 = Nussinov(rnas50[i])
        F100, score100 = Nussinov(rnas100[i])
        sums = [sum(x) for x in zip(sums, [score10, score50, score100])]
    
    ## PLOT AVERAGE SCORE VS SEQ LENGTH

    averageScores = [sums[0]/numSamples, sums[1]/numSamples, sums[2]/numSamples]
    plt.plot([10, 50, 100], averageScores)
    plt.axis([0, 110, 0, max(averageScores)+10])
    plt.ylabel("Average score")
    plt.xlabel("Length of RNA sequence")
    plt.show()


# part D
def variousGCcontent():

    ## GENERATE RANDOM RNAS BY UNIFORM DISTRIBUTION

    numSamples = 1000
    rnas100 = []
    gcContent = []
    for i in range(numSamples):
        gc = float(i % 10) / 10.0
        rnas100.append(GenerateRandomRna(100, gc))
        gcContent.append(gc)

    ## RUN NUSSINOV ALGORITHM

    scores = []   # part D 
    for i in range(numSamples):
        F, s = Nussinov(rnas100[i])
        scores.append(s)   # part D
    
    ## PLOT AVERAGE SCORE VS SEQ LENGTH

    plt.scatter(x=gcContent, y=scores)
    plt.ylabel("Average score")
    plt.xlabel("GC content")
    plt.show()


##### MAIN #####

def main():

    variousGCcontent()



if __name__ == "__main__":
    main()