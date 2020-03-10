#!/usr/bin/env python

import sys

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3

# ORIGINAL CODE for PART A with one section filled in
def seqalignDP_PartA(seq1,seq2,subst_matrix,gap_pen):
	"""return the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
	   Note: gap_pen should be positive (it is subtracted)
	"""
	F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
	TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

	#	S2 --> 0,..,j,..,len(S2)+1
	# S1
	# |
	# V
	# 0,..,i,..,len(S1)+1

	# initialize dynamic programming table for Needleman-Wunsch alignment (Durbin p.20)
	# initializes first row and first column

	# rows
   	for i in range(1,len(seq1)+1):
   		F[i][0] = 0 - i*gap_pen
		TB[i][0] = PTR_GAP2 # indicates a gap in seq2

	# cols
	for j in range(1,len(seq2)+1):
   		F[0][j] = 0 - j*gap_pen
		TB[0][j] = PTR_GAP1 # indicates a gap in seq1

	# YOUR CODE HERE
	# Fill in the dynamic programming tables F and TB, starting at [1][1]
	# Hints: The first row and first column of the table F[i][0] and F[0][j] are dummies
	#        (see for illustration Durbin p.21, Figure 2.5, but be careful what you
	#         think of as rows and what you think of as columns)
	#        Hence, the bases corresponding to F[i][j] are actually seq1[i-1] and seq2[j-1].
	#        Use the dictionary base_idx to convert from the character to an index to
	#         look up entries of the substitution matrix.
	#        To get started, you can complete and run the algorithm filling in only F,
	#         and then figure out how to do TB.

	for i in range (1, len(seq1)+1):

		for j in range (1, len(seq2)+1):

			# calculate base pair score
			bpScore = subst_matrix[base_idx.get(seq1[i-1])][base_idx.get(seq2[j-1])]

			# find max of all previous scores and set traceback
			score1 = F[i-1][j] - gap_pen
			score2 = F[i][j-1] - gap_pen
			score3 = F[i-1][j-1] + bpScore
			F[i][j] = score3
			TB[i][j] = PTR_BASE
			if score1 > F[i][j]:
				F[i][j] = score1
				TB[i][j] = PTR_GAP2
			if score2 > F[i][j]:
				F[i][j] = score2
				TB[i][j] = PTR_GAP1

	return F[len(seq1)][len(seq2)], F, TB


# PART B, used a min distance type of score, such that all scores > 0
# a seq aligned with itself should have a zero score,
# greater mismatch confers a higher score
def seqalignDP_PartB(seq1,seq2,subst_matrix,gap_pen):
	"""return the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
	   Note: gap_pen should be positive (it is subtracted)
	"""
	F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
	TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

	#	S2 --> 0,..,j,..,len(S2)+1
	# S1
	# |
	# V
	# 0,..,i,..,len(S1)+1

	# initialize dynamic programming table for Needleman-Wunsch alignment (Durbin p.20)
	# initializes first row and first column

	# rows
   	for i in range(1,len(seq1)+1):
   		F[i][0] = 0 + i*gap_pen
		TB[i][0] = PTR_GAP2 # indicates a gap in seq2

	# cols
	for j in range(1,len(seq2)+1):
   		F[0][j] = 0 + j*gap_pen
		TB[0][j] = PTR_GAP1 # indicates a gap in seq1

	# fill out table F
	for i in range (1, len(seq1)+1):

		for j in range (1, len(seq2)+1):

			# calculate base pair score
			bpScore = subst_matrix[base_idx.get(seq1[i-1])][base_idx.get(seq2[j-1])]
			# if seq1[i-1] == "T" and seq2[j-1] == "G":
			# 	print "i = ", i, ", S1[i] = ", seq1[i-1]
			# 	print "j = ", j, ", S2[j] = ", seq2[j-1]
			# 	print bpScore
			# 	print "--------"

			# # find max of all previous scores and set traceback
			score1 = F[i-1][j] + gap_pen
			score2 = F[i][j-1] + gap_pen
			score3 = F[i-1][j-1] + bpScore**2
			F[i][j] = score3
			TB[i][j] = PTR_BASE
			if score1 < F[i][j]:
				F[i][j] = score1
				TB[i][j] = PTR_GAP2
			if score2 < F[i][j]:
				F[i][j] = score2
				TB[i][j] = PTR_GAP1

	return F[len(seq1)][len(seq2)], F, TB


def traceback(seq1,seq2,TB):
	s1 = ""
	s2 = ""

	i = len(seq1)
	j = len(seq2)

	while TB[i][j] != PTR_NONE:
		if TB[i][j] == PTR_BASE:
			s1 = seq1[i-1] + s1
			s2 = seq2[j-1] + s2
			i=i-1
			j=j-1
		elif TB[i][j] == PTR_GAP1:
			s1 = '-' + s1
			s2 = seq2[j-1] + s2
			j=j-1
	   	elif TB[i][j] == PTR_GAP2:
			s1 = seq1[i-1] + s1
			s2 = '-' + s2
			i=i-1
		else: assert False

	return s1,s2

def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):		# if string starts with > character, skip line
            continue
        seq.append(line.rstrip())		# else, append line to sequence and remove trailing characters, default is spaces
	
    return "".join(seq)					# converts array seq to string and returns

def checkIfEqual(s1, s2):

	isOkay = True
	mismatches = []

	if (len(s1) != len(s2)):
		isOkay = False
		mismatches.append(0)
		return isOkay, mismatches
	
	for i in range(0,len(s1)):
		if s1[i] == s2[i]:
			continue
		elif s1[i] == "-" or s2[i] == "-":
			continue
		else:
			isOkay = False
			mismatches.append(i)
	
	return isOkay, mismatches


## my scoring matrix, taken from the lecture, PART A
S_partA = [
	# A  G   C   T
	[ 1,  0, -1, -1], # A
	[ 0,  1, -1, -1], # G
	[-1, -1,  1,  0], # C
	[-1, -1,  0,  1]  # T
]

gap_pen_partA = -1

## my scoring matrix, PART B
S_partB = [
	# A  G   C   T
	[ 0,  2,  3,  3], # A
	[ 2,  0,  3,  3], # G
	[ 3,  3,  0,  2], # C
	[ 3,  3,  2,  0]  # T
]

gap_pen_partB = 1

## their scoring matrix
## found gap penalty was too high, led to incorrect traceback
# S = [
# 	# A  G   C   T
# 	[3, -1, -2, -2], # A
# 	[-1, 3, -2, -2], # G
# 	[-2, -2, 3, -1], # C
# 	[-2, -2, -1, 3]  # T
# 	]
# gap_pen = 4

def main():

    # parse commandline
	# if len(sys.argv) < 3:
	# 	print "you must call program as: python ps1-seqalign.py <FASTA 1> <FASTA 2>"
	# 	sys.exit(1)

	# file1 = sys.argv[1]
	# file2 = sys.argv[2]
	file1 = "human_HoxA13.fa"
	file2 = "mouse_HoxA13.fa"

	seq1 = readSeq(file1)
	seq2 = readSeq(file2)

	## PART A
	# score, F, TB = seqalignDP_PartA(seq1,seq2,S_partA,gap_pen_partA)

	## PART B
	score, F, TB = seqalignDP_PartB(seq1,seq2,S_partB,gap_pen_partB)

	print >> sys.stderr, "Score = ", score

	s1, s2 = traceback(seq1,seq2,TB)
	# print s1
	# print s2

	isOkay, mismatches = checkIfEqual(s1, s2)

	if isOkay:
		print "Alignment successful!"
	else:
		print "Alignment failed - ", len(mismatches), " mismatches detected"
		m = mismatches[0] + 1
		print "First mismatch at ", m
		subs1 = []
		subs2 = []
		for i in range(-2,3):
			subs1.append(s1[m-1+i])
			subs2.append(s2[m-1+i])
		print subs1
		print subs2
		print "--------"
		print "F[i] = "
		print F[m-1][m-1], " | ", F[m-1][m], " | " , F[m-1][m+1]
		print "-" * 15
		print F[m][m-1], " | ", F[m][m], " | " , F[m][m+1]
		print "-" * 15
		print F[m+1][m-1], " | ", F[m+1][m], " | " , F[m+1][m+1]
		

if __name__ == "__main__":
	main()
