# Problem 2, parts B and C
# My understanding of this problem is that we are supposed to generate
# 10,000 alignments (each alignment is 4 sequences) of length 10 
# using the conditional probability  distribution of conserved and 
# non-conserved regions given in part A.

# Score is the number of unique pairs that share the same symbol

#   Score           P(S|N)      P(S|C)
#  -------------------------------------
#   0               0.1         0.05
#   1               0.35        0.15
#   2               0.25        0.2
#   3               0.2         0.3
#   6               0.1         0.3

import numpy as np

base_char = { 0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T' }
base_idx = {'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
score_idx = { 0 : 0, 1 : 1, 2 : 2, 3 : 3, 6 : 4 }
state_idx = { 'N' : 0, 'C' : 1 }

# the probability distribution of scores, or the 
# Probability Mass Function (PMF)
condprob = [
    # score -> 0, 1, 2, 3, 6
    [0.1, 0.35, 0.25, 0.2, 0.1],    # N
    [0.05, 0.15, 0.2, 0.3, 0.3]     # C
]

N = 10  # length of alignments
M = 4   # number of sequences per alignment
S = 5   # number of possible score outcomes
B = 4   # number of bases

# Example of an alignment:

    # ACGACGACTA
    # CAGACGCTGA
    # TTCCTCTGAT
    # AGATGTGACT

# In this example, N = 10, M = 4. B and S are not affected by definition of an alignment

def getRandomBase(skip=[]):
    
    choices = []                # allowed base choices
    for i in range(B):
        if base_char[i] not in skip:
            choices.append(base_char[i])
    
    if len(choices) == 0:       # if all bases are in skip, return 'X'
        return 'X'
    elif len(choices) == 1:     # if only one choice, just return that one
        return choices[0]

    u = np.random.uniform()     # random sample from uniform distribution
    udiv = 1.0/len(choices)

    for i in range(1, len(choices)+1):
        if i <= len(choices) and u < i*udiv:
            return choices[i-1]
    
    # if didn't return because udiv*(len(choices) + 1) < u < 1.0,
    # return last choice
    return choices[len(choices) - 1]


def sampleProbabilityDistribution(state):

    if state not in state_idx:
        print state, " is not a valid state"
        return []

    iState = state_idx[state]

    # find the cumulative probability distribution
    cdf = []
    for i in range(S):
        if i == 0:
            cum_prob = condprob[iState][i]
        else:
            cum_prob = condprob[iState][i] + cdf[i-1]
        cdf.append(cum_prob)
    
    # generate a uniform distribution
    U = np.random.uniform(0,1,N)

    # generate samples from cdf using the uniform samples
    sampleScores = []
    for u in U:
        if u < cdf[0]:
            sampleScores.append(0)
        elif u < cdf[1]:
            sampleScores.append(1)
        elif u < cdf[2]:
            sampleScores.append(2)
        elif u < cdf[3]:
            sampleScores.append(3)
        else:
            sampleScores.append(6)
    
    return sampleScores


def generateAlignment(scores):

    if len(scores) != N:
        print "Size of scores is not equal to ", N
        return []

    seqs = []
    for i in range(M):
        seqs.append([])

    for s in scores:

        if s == 6:                    # Easy case: all bases are the same

            b = getRandomBase()
            for i in range(M):
                seqs[i].append(b)


        elif s == 0:                  # Easy case: no pairs, each base is different

            bases = []
            for i in range(M):
                b = getRandomBase(bases)
                bases.append(b)
                seqs[i].append(b)
        
        else:

            # let numSameBases be the reoccuring number of bases belonging to pairs in an 
            # alignment column
            #   example, AATG -> 2, AATT -> 2, AAAT -> 3
            # let numDiffPairs be the number of unique BP pairs in an alignment column
            #   example, AATG -> 1, AATT -> 2, AAAT -> 1
            # let numDiffBases be the number of different BPs in an alignment column
            #   example, AATG -> 3, AATT -> 2, AAAT -> 2, 

            if s == 1:
                numSameBases = 2
                numDiffPairs = 1
                numDiffBases = 3
            elif s == 2:
                numSameBases = 2
                numDiffPairs = 2
                numDiffBases = 2
            elif s == 3:
                numSameBases = 3
                numDiffPairs = 1
                numDiffBases = 2
            
            i = 0
            bases = ['X' for bp in range(M)]

            while i < numDiffBases:

                b = getRandomBase(bases)

                if i < numDiffPairs: 
                    j = 0          
                    while j < numSameBases:
                        pos = np.random.random_integers(0,M-1)
                        if bases[pos] == 'X':
                            bases[pos] = b
                            j += 1   
                else:
                    for k in range(M):
                        if bases[k] == 'X':
                            bases[k] = b  
                            break              
                
                i += 1
            
            for sq in range(M):
                seqs[sq].append(bases[sq])



    alignment = []
    for sq in range(B):
        alignment.append("".join(seqs[sq]))

    return alignment


def calculateConditionalProbabilities(scores):

    if len(scores) != N:
        print "Sequences are the incorrect length"
        return 0.0, 0.0
        
    pN, pC = 1.0, 1.0

    for i in range(N):
        pN *= condprob[state_idx['N']][score_idx[scores[i]]]
        pC *= condprob[state_idx['C']][score_idx[scores[i]]]
    
    return pN, pC


if __name__ == "__main__":

    numSamples = 10000
    state = 'N'
    counts = [0, 0] # [num where N > C, num where C > N]

    for i in range(numSamples):

        # generate random scores given probability distribution for N and C
        scores = sampleProbabilityDistribution(state)

        # not necessary to find conditional probabilities,
        # but did it anyway
        # seqs = generateAlignment(scores)

        pN, pC = calculateConditionalProbabilities(scores)

        if (pN > pC):
            counts[state_idx['N']] += 1
        elif (pN < pC):
            counts[state_idx['C']] += 1
        
    # print results
    print "State selected was ", state
    print "% of samples where P(S|N) > P(S|C) = ", (float(counts[state_idx['N']]) / numSamples) * 100
    print "% of samples where P(S|C) > P(S|N) = ", (float(counts[state_idx['C']]) / numSamples) * 100


        