## Using the same HMM as PS1/viterby.py

## Posterior decoding calculates the total probability that emission xi came from
## state k across all paths. It gives us the path containing the most likely state
## at any point in time (pi^ = argmax(k)P{pi_i = k | x})

## Different from Viterbi, which finds the optimal path given a series of emisions, X
## (finds max scoring path over all paths, pi* = argmax(pi)P{x | pi})

## Reference 
## https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm#Python_example

import sys
from math import log
import numpy as np

###############################################################################
# HMM PARAMETERS
# Conventions: + refers to High-GC, and - refers to Low-GC. When indexing
#  state_idx, 0 is + and 1 is -.
###############################################################################

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
state_idx = { '+' : 0, '-' : 1 }

# initial distribution over state_idx, i.e. probability of starting in state k
init_dist = [0.5,0.5]
final_dist = [0.5,0.5]

# transition probabilities
tr = [
    #  to+   to-
    [ 0.99, 0.01 ], # from+
    [ 0.01, 0.99 ]  # from-
]


# emission probabilities
em = [
    #    A     G     C     T
    [ 0.20, 0.30, 0.30, 0.20], # +
    [ 0.30, 0.20, 0.20, 0.30]  # -
]


def posteriorDecode(sequence):

    fmatrix = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]
    bmatrix = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]    

    # Forward part of the algorithm
    for i in range(0, len(sequence)):    	
        iBase = base_idx[sequence[i]]        
        for st in state_idx:
            iState = state_idx[st]
            if i == 0:
                fsumPrev = log(init_dist[iState])
            else:
                fsumPrev = sum(fmatrix[state_idx[k]][i-1] + log(tr[state_idx[k]][iState]) for k in state_idx)

            fmatrix[iState][i] = log(em[iState][iBase]) + fsumPrev

    # calculate the probability of all emissions
    log_p_fwd = sum(fmatrix[state_idx[k]][len(sequence)-1] + log(final_dist[state_idx[k]]) for k in state_idx)
    # convert p_fwd from log to 10^x
    p_fwd = 10**log_p_fwd

    # Backward part of the algorithm
    for i in range(len(sequence)-1, -1, -1):
        iBase = base_idx[sequence[i]]
        for st in state_idx:
            iState = state_idx[st]
            if i == len(sequence)-1:
                bmatrix[iState][i] = final_dist[iState]
            else:
                iBase = base_idx[sequence[i+1]] 
                bmatrix[iState][i] = sum(em[state_idx[l]][iBase]*tr[state_idx[l]][iState]*bmatrix[state_idx[l]][iBase] for l in state_idx)

    # iBase = base_idx[sequence[0]]
    # p_bkw = sum(init_dist[l] * em[l][iBase] * bmatrix[state_idx[l]][iBase] for l in state_idx)

    # Merging the two parts
    posterior = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]
    if p_fwd == 0:
        print "P_fwd is zero! ", p_fwd
        return False, fmatrix, bmatrix, posterior
    for i in range(0, len(sequence)):
        for st in state_idx:
            iState = state_idx[st]
            iBase = base_idx[sequence[i]]
            posterior[iState][i] = fmatrix[iState][i] * bmatrix[iState][i] / p_fwd

    return True, fmatrix, bmatrix, posterior

def main():

    # parse command-line arguments
    # if len(sys.argv) < 2:
    #     print "you must call program as: ./posterior_decoding.py <datafile>"
    #     sys.exit(1)
    # datafile = sys.argv[1]
    datafile = "PS1/hmmgen"

    # read sequences
    print "reading sequence"
    f = open(datafile)
    X = f.readline()
    # refanno = f.readline()
    f.close()

    if X[len(X)-1] == '\n': X=X[0:len(X)-1]
    # if refanno[len(refanno)-1] == '\n': refanno=refanno[0:len(refanno)-1]

    X=list(X)
    # refanno=list(refanno)
    # for i in xrange(len(X)):
    #     X[i] = base_idx[X[i]]
    # for i in xrange(len(refanno)):
    #     refanno[i] = state_idx[refanno[i]]

    # perform posterior decoding
    print "performing posterior decoding"
    success, fmatrix, bmatrix, posterior = posteriorDecode(X)

    if not success:
        print "PD failed, exiting now"
        exit(1)

    # find the path containing the most likely state at each position
    print "finding path"
    N = len(posterior[0])
    pi = []
    for i in range(0, N):
        # find the state with the max probability and add to pi
        maxState = '+'   
        for st in state_idx:
            if st != maxState and posterior[iState][i] > posterior[state_idx[maxState]][i]:
                maxState = st
        pi.append(maxState)
    
    print "-" * 15
    print "The path containing the most likely state at each position is"
    print pi

if __name__ == "__main__":
    main()