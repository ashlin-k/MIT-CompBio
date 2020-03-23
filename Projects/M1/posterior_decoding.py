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
    [ 0.80, 0.20 ], # from+
    [ 0.20, 0.80 ]  # from-
]


# emission probabilities
em = [
    #    A     G     C     T
    [ 0.20, 0.30, 0.30, 0.20], # +
    [ 0.30, 0.20, 0.20, 0.30]  # -
]

###############################################################################
# POSTERIOR DECODING ALGORITHM 
###############################################################################

def posteriorDecode(sequence):

    # init matrices

    fmatrix = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]
    bmatrix = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]  

    fmatrix[0][0] = log(init_dist[0])  
    fmatrix[1][0] = log(init_dist[1])  
    bmatrix[0][len(sequence)-1] = log(final_dist[0])  
    bmatrix[1][len(sequence)-1] = log(final_dist[1])  

    # Forward part of the algorithm
    for i in range(1, len(sequence)):    	
        iBase = base_idx[sequence[i]]        
        for st in state_idx:
            iState = state_idx[st]
            # print "fmatrix[st = 0] = ", fmatrix[0][i-1], "log(tr) = ", log(tr[0][iState]), \
            #     "fmatrix[st = 1] = ", fmatrix[1][i-1], "log(tr) = ", log(tr[1][iState])                
            fsumPrev = sum(fmatrix[state_idx[k]][i-1] + log(tr[state_idx[k]][iState]) for k in state_idx)
            fmatrix[iState][i] = log(em[iState][iBase]) + fsumPrev
            # print "sum = ", fmatrix[iState][i]
            # print "-----------"

    # calculate the probability of all emissions
    log_p_fwd = sum(fmatrix[state_idx[k]][len(sequence)-1] + log(final_dist[state_idx[k]]) for k in state_idx)
    # convert p_fwd from log to 10^x
    p_fwd = 10**log_p_fwd

    # Backward part of the algorithm
    for i in range(len(sequence)-2, -1, -1):
        iBase = base_idx[sequence[i]]
        for st in state_idx:
            iState = state_idx[st]
            iBase = base_idx[sequence[i+1]] 
            bmatrix[iState][i] = sum(log(em[state_idx[l]][iBase]) + \
                log(tr[state_idx[l]][iState]) + bmatrix[state_idx[l]][iBase] for l in state_idx)

    # iBase = base_idx[sequence[0]]
    # p_bkw = sum(init_dist[l] * em[l][iBase] * bmatrix[state_idx[l]][iBase] for l in state_idx)

    # Merging the two parts
    posterior = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]
    # if p_fwd == 0:
    #     print "P_fwd is zero! ", p_fwd
    #     return False, fmatrix, bmatrix, posterior
    for i in range(0, len(sequence)):
        for st in state_idx:
            iState = state_idx[st]
            iBase = base_idx[sequence[i]]
            # posterior[iState][i] = fmatrix[iState][i] * bmatrix[iState][i] / log_p_fwd
            posterior[iState][i] = fmatrix[iState][i] + bmatrix[iState][i] - log_p_fwd

    return True, fmatrix, bmatrix, posterior


def getPath(fmatrix, bmatrix, posterior):

    N = len(posterior[0])
    pi = []
    for i in range(0, N):
        # find the state with the max probability and add to pi
        maxState = '+'   
        for st in state_idx:
            iState = state_idx[st]
            if st != maxState and posterior[iState][i] > posterior[state_idx[maxState]][i]:     # comparing log values
                maxState = st
        pi.append(state_idx[maxState])
    
    return pi

###############################################################################
# ANNOTATION BENCHMARKING
###############################################################################

def basecomp(X,anno):
    counts = [[0]*4,[0]*4]
    for i in xrange(len(X)):
        counts[anno[i]][X[i]] += 1
    sum0 = sum(counts[0])
    sum1 = sum(counts[1])

    basecomp0, basecomp1 = [0 for i in range(len(base_idx))], [0 for i in range(len(base_idx))]
    if sum0 != 0:
        basecomp0 = [float(x)/sum0 for x in counts[0]]
    if sum1 != 0:
        basecomp1 = [float(x)/sum1 for x in counts[1]]
    
    return [basecomp0, basecomp1]

def region_lengths(anno):
    lengths = [[],[]]
    curlen=1
    for i in xrange(1,len(anno)):
        if anno[i] == anno[i-1]:
            curlen += 1
        else:
            lengths[anno[i-1]].append(curlen)
            curlen=1
    lengths[anno[len(anno)-1]].append(curlen)
    return lengths

def anno_accuracy(refanno,testanno):
    correct = 0
    assert len(refanno) == len(testanno)
    for i in xrange(1,len(refanno)):
        if refanno[i] == testanno[i]:
            correct += 1
    return float(correct)/len(refanno)

def print_basecomp(b):
    print "A=%.2f%% G=%.2f%% C=%.2f%% T=%.2f%%" % (100*b[0],100*b[1],100*b[2],100*b[3])

def print_annostats(X,anno,filename):
    lengths = region_lengths(anno)
    basecomps = basecomp(X,anno)

    highGClength, lowGClength = 0, 0
    if len(lengths[0]) > 0:
        highGClength = sum(lengths[0])/len(lengths[0])
    if len(lengths[1]) > 0:
        lowGClength = sum(lengths[1])/len(lengths[1])

    print "High-GC mean region length: ", highGClength
    print "High-GC base composition:",
    print_basecomp(basecomps[0])
    print "Low-GC mean region length: ", lowGClength
    print "Low-GC base composition:",
    print_basecomp(basecomps[1])


###############################################################################
# MAIN
###############################################################################

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
    refanno = f.readline()
    f.close()

    if X[len(X)-1] == '\n': X=X[0:len(X)-1]
    if refanno[len(refanno)-1] == '\n': refanno=refanno[0:len(refanno)-1]

    Xlist=list(X)
    refanno=list(refanno)
    for i in xrange(len(Xlist)):
        Xlist[i] = base_idx[X[i]]
    for i in xrange(len(refanno)):
        refanno[i] = state_idx[refanno[i]]

    # perform posterior decoding
    # all matrices encode the log value of probabilities
    print "performing posterior decoding"
    success, fmatrix, bmatrix, posterior = posteriorDecode(X)

    if not success:
        print "PD failed, exiting now"
        exit(1)

    # find the path containing the most likely state at each position
    print "finding path"
    pi = getPath(fmatrix, bmatrix, posterior)
    
    print ""

    print "Authoritative annotation statistics"
    print "-----------------------------------"
    print_annostats(Xlist,refanno,datafile+"_authoritative")
    print ""

    print "PD annotation statistics"
    print "-----------------------------"
    print_annostats(Xlist,pi,datafile+"_PD")
    print ""
    
    print "Accuracy: %.2f%%" % (100*anno_accuracy(refanno,pi))

if __name__ == "__main__":
    main()