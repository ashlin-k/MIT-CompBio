# Problem 2, parts B and C
# My understanding of this problem is that we are supposed to generate
# 10,000 sequences of length 10 using the conditional probability 
# distribution of conserved and non-conserved regions given in part A.
# Score is the number of unique pairs that share the same symbol

#   Score           P(S|N)      P(S|C)
#  -------------------------------------
#   0               0.1         0.05
#   1               0.35        0.15
#   2               0.25        0.2
#   3               0.2         0.3
#   6               0.1         0.3

base_idx = { 0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T' }
score_idx = { 0 : 0, 1 : 1, 2 : 2, 3 : 3, 6 : 4 }
state_idx = { 'N' : 0, 'C' : 1 }

condprob = [
    # score -> 0, 1, 2, 3, 6
    [0.1, 0.35, 0.25, 0.2, 0.1],    # N
    [0.05, 0.15, 0.2, 0.3, 0.3]     # C
]

N = 10  # length of alignments

def generateSequence(state):

    for i in range(N):

        