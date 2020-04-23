## PS3 problem 3: RNA folding, Nussinov algorithm

import numpy as np

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'U' : 3}

# base pair scoring matrix
S = [
    #   A   G   C   U
    [   0,  0,  0,  -1],    # A
    [   0,  0,  -1, -1],    # G
    [   0,  -1, 0,  0 ],    # C
    [   -1, -1, 0,  0 ]     # U
]

# store scores in F
# F = [
#     1 2 . . . N
#   [. . . . . . .],    1
#   [. . . . . . .],    2
#   . . .
#   [. . . . . . .]     N
# ]


def calculateScore(i, j, seq, F):

    Si = seq[i]
    Sj = seq[j]
    score1 = F[i+1][j-1] + S[base_idx[Si]][base_idx[Sj]]
    score2 = [1+ F[i][k] + F[k+1][j] for k in range(i, j)]
    return max(score1, max(score2))


def Nussinov(seq):

    N = len(seq)    
    F = [ [0.0 for i in range(N+1)] for j in range(N+1)]

    for l in range(2, N+1):
        for i in range(1, N+1-l):
            j = i + l - 1
            F[i][j] = calculateScore(i, j, seq, F)
    
    return F


def GenerateRandomRna(n):

    rna = ""
    for i in range(n):
        r = np.random.uniform()
        if r < 0.25:
            rna += 'A'
        elif r < 0.5:
            rna += 'G'
        elif r < 0.75:
            rna += 'C'
        else:
            rna += 'U'
    
    return rna

def main():

    # rnas = []
    # for i in range(1000):
    #     rnas.append(GenerateRandomRna(100))
    
    rna = GenerateRandomRna(10)

    F = Nussinov(rna)

    for i in range(len(F)):
        print F[i]
    print ""

    


if __name__ == "__main__":
    main()