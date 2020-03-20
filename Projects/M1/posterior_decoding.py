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
from util import plothist

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


def fwd_bkw(sequence):
    """Forward–backward algorithm."""

    fmatrix = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]
    bmatrix = [[0 for col in range(len(sequence))] for row in range(len(state_idx))]

    # Forward part of the algorithm
    for i in range(len(sequence)):    	
    	iBase = base_idx[sequence[i]]        
        for st in state_idx:
        	iState = state_idx[st]
            if i == 0:
                fsumPrev = init_dist[iState]
            else:
                fsumPrev = sum(fmatrix[k][i-1]*tr[k][iState] for k in state_idx)

            fmatrix[iState][i] = em[iState][iBase] * fsumPrev


    p_fwd = sum(fmatrix[k][len(sequence)-1] * final_dist[state_idx[k]] for k in state_idx) # ??? we don't know where it will end

    # Backward part of the algorithm






    bkw = []
    b_prev = {}
    for i, observation_i_plus in enumerate(reversed(observations[1:]+(None,))):
        b_curr = {}
        for st in state_idx:
            if i == 0:
                # base case for backward part
                b_curr[st] = tr[st][end_st]
            else:
                b_curr[st] = sum(tr[st][l] * em[l][observation_i_plus] * b_prev[l] for l in state_idx)

        bkw.insert(0,b_curr)
        b_prev = b_curr

    p_bkw = sum(init_dist[l] * em[l][observations[0]] * b_curr[l] for l in state_idx)

    # Merging the two parts
    posterior = []
    for i in range(len(observations)):
        posterior.append({st: fwd[i][st] * bkw[i][st] / p_fwd for st in state_idx})

    assert p_fwd == p_bkw
    return fwd, bkw, posterior