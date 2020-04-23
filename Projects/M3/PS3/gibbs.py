#!/usr/bin/env python
import sys
import string
import random
import math
import numpy as np

# Reference:
# https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gibbs.pdf

#### INSTRUCTIONS FOR USE:
# call program as follows: ./gibbs.py <Motif Length> <Data File>
# make sure the gibbs.py is marked as executable: chmod +x gibbs.py

alphabet = ['A', 'G', 'C', 'T']
base_idx = {'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3}

# initialize base frequencies
B = [0.25 for a in range(len(alphabet))]   

#### GibbsSampler:
#### 	INPUTS:	S - list of sequences
####		L - length of motif
####	OUTPUT:	PWM - 4xL list with frequencies of each base at each position
####                  Order of bases should be consistent with alphabet variable
####              |     0   1   2   ...     L-1
####            ----------------------------------
####            A | [[                          ],
####            G |  [                          ],
####            C |  [                          ],
####            T |  [                          ]]

def GibbsSampler(S,L):

    # initialize the Position Weight Matrix (PWM)
    # set all probabilities equal
    PWM = []
    for i in range(len(alphabet)):
        PWM.append([1.0/len(alphabet)]*L) 

    ######### ADD YOUR CODE HERE ######

    # initialize random starting points for all seqs in S
    A = [random.randint(0, len(sq)-L) for sq in S]
    prevA = None

    while A != prevA:

        prevA = A

        for i in range(len(S)):

            # create a profile of motifs L long, for all seqs except S[i]
            # we'll be checking S[i] against all motifs in motifProfiles
            motifProfiles = [sq[j : j + L] for q, (sq,j) in enumerate(zip(S, A)) if q != i]

            # recalculate PWM model
            PWM = UpdateModel(motifProfiles)

            # find highest scoring start index in S[i]
            iBest = None
            bestScore = -1
            for j in range(0, len(S[i]) - L + 1):
                kmer = S[i][j : j+L]
                score = Score(PWM, kmer)
                if score > bestScore or iBest == None:
                    iBest = j
                    bestScore = score
            
            A[i] = iBest
    
    return PWM

    ######### END OF YOUR CODE HERE #####
    

###### YOUR OWN FUNCTIONS HERE
# optional -- feel free to add your own functions if you want to

def Score(M, seq):
    
    if len(M) != len(alphabet):
        print "Incorrect dimensions of M"
        return -1

    if len(M[0]) != len(seq):
        print "Num of columns of M is not equal to the length of the sequence"
        return -1

    L = len(seq)
    score = 0.0
    for i in range(L):
        iBase = base_idx[seq[i]]
        if M[iBase][i] != 0:
            score += math.log(M[iBase][i] / B[iBase])
        else:   # if M is 0, use smallest possible number > 0
            score += math.log(np.nextafter(0, 1))
    
    return score


def UpdateModel(profiles):
    
    L = len(profiles[0])
    N = float(len(profiles))

    PWM = []
    for i in range(len(alphabet)):
        PWM.append([0]*L) 

    # first get counts
    for i in range(len(profiles)):
        motif = profiles[i]
        for j in range(L):
            PWM[base_idx[motif[j]]][j] += 1 
    
    # second express as a probability
    for j in range(L):
        for a in range(len(alphabet)):
            PWM[a][j] /= N
    
    return PWM


def GetMotif(PWM):

    motif = []
    L = len(PWM[0])
    
    for j in range(L):
        bestP = 0.0
        bestBase = 'X'
        for base in base_idx:
            iBase = base_idx[base]
            prob = PWM[iBase][j]
            if prob > bestP or bestBase == 'X':
                bestP = prob
                bestBase = base
        motif.append(bestBase)
    
    return "".join(motif)

def PrintResults(P, motif):

    L = len(P[0])

    # print PWM
    print "    ", 
    for i in range(L):
        print "%-5d " % (i+1),
    print ""
	
    for j in range(len(alphabet)):
        print " %s " % alphabet[j], 
        for i in range(L):
            print " %5.3f" % P[j][i],
        print ""
    
    # print motif
    print "-" * (4+7*L)
    print "    ",
    for i in range(L):
        print "%s     " % motif[i],
    print ""
    print ""


###### END OF YOUR FUNCTIONS

def readdata(file):
    data = []
    for line in open(file,'r'):
        data.append(line[0:-1])
    return data


def main():

    L = int(sys.argv[1])
    datafile = sys.argv[2]
    # L = 3
    # datafile = "gibbs_simple_test.txt"
    S = readdata(datafile)
	
    # run gibbs sampling multiple times to find the mostly likely motif
    N = 10  # 100
    motifDict = {}
    for i in range(N):
        P = GibbsSampler(S,L)
        motif = GetMotif(P)
        motifDict.setdefault(motif, []).append(P)
        PrintResults(P, motif)
    
    ### EXTRA: FIND THE MOST COMMON PWM PROGRAMMATICALLY ###

    # from iterations, find the average P matrix
    mostCommonMotif = ""
    count = -1
    for motif in motifDict:
        if len(motifDict[motif]) > count:
            mostCommonMotif = motif
            count = len(motifDict[motif])
    
    # check we found at least one motif
    if count == -1 or mostCommonMotif == "":
        print "No motifs found"
        exit(1)
    
    # find the most common P for the mostCommonMotif
    Pcounts, Phashes = {}, {}
    for i in range(count):
        Pi = motifDict[mostCommonMotif][i]
        Phash = hash(str(Pi))
        if Phash not in Phashes:
            Phashes[Phash] = Pi
            Pcounts[Phash] = 1
        else:
            Pcounts[Phash] += 1
    
    # choose the P matrix with the highest count
    count = -1
    bestHash = 0
    for h in Pcounts:
        if Pcounts[h] > count:
            count = Pcounts[h]
            bestHash = h
    
    # print "#" * (4+7*L)
    # print ""
    print "The most common PWM is:"
    print ""
    PrintResults(Phashes[bestHash], mostCommonMotif)
	


if __name__ == "__main__":
    main()