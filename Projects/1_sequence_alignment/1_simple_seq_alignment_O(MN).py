#	For a sequence S1 of length M, and a sequence S2 of length N,
#	this is a simple O(M*N), in both time and space, 
#	algorithm for in exact DNA sequence alignment.
#	This is a dynamic programming algorithm which scores base pair 
# 	alignment and gaps in order to find the lowest cost / min change
#	alignment.

#	Based of of Lecture 2 in the MIT CompBio lecture series
#	https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-047-computational-biology-fall-2015/lectures_slides/MIT6_047F15_Lecture02.pdf


## def scorePair	This function scores the cost of aligning two base pairs 
##					according to the below matrix.  direct match is 1, a 
##					changed from a purine to another purine (A,G) is 0,
##					a change pyrimidine to a pyrimidine (T, C) is also 0,
##					and a change from a purine to a pyrimidine and vice versa
##					is -1. Introducing a gap is -1.

##		A 	G 	T 	C
##	A 	1	0	-1	-1
##	G 	0 	1 	-1 	-1
##	T 	-1	-1	1 	0
##	C 	-1 	-1 	0 	1
## 	gap = -1

## arg p1			char, first base pair
##					'A', 'G', 'T', 'C'
## arg p2			char, second base pair
##					'A', 'G', 'T', 'C'

def scorePair(p1, p2): 			## O(1)

	if (p1 == p2):
		return 1

	elif (p1 == 'A' or p1 == 'G') and (p2 == 'A' or p2 == 'G'):
		return 0

	elif (p1 == 'T' or p1 == 'C') and (p2 == 'T' or p2 == 'C'):
		return 0

	else:
		return -1


## def scoreCell	For the 2D data matrix which stores the scores of all 
##					combos of sequence alignments, this function calculates 
##					the score of cell (i, j) based on the values of 3 previous
##					cells:
##		score(i,j) = max {
##			score(j-1, i) - gap,
##			score(j, i-1) - gap,
##			score(j-1,i-1) + scorePair(S1[i], S2[j])
##		}

## arg i 			int, index of BP in S1; indicates the column in datamatrix
## arg j 			int, index of BP in S2; indicates the row in datamatrix	

def scoreCell(i, j):

	global M, N, datamatrix, S1, S2

	if (i == 0 and j == 0):		
		return 0;

	gap = -1

	score, score2, score3 = -10, -10, -10

	if (i > 0 and j > 0):
		score = datamatrix[j-1][i-1] + scorePair(S1[i-1], S2[j-1])

	if (i > 0):
		score2 = datamatrix[j][i-1] + gap
		if (score2 > score):
			score = score2

	if (j > 0):
		score3 = datamatrix[j-1][i] + gap
		if (score3 > score):
			score = score3

	datamatrix[j][i] = score

##	def findOptimalSequence		After datamatrix has been filled out,
##								this function goes traces backward 
##								(from bottom-right to top-left), the highest
##								scoring sequence, and then creates the new 
##								S1 and S2 sequences, inserting gaps where
##								horizontal and vertical transitions occur.
##		horizontal transition -> i changes, j is the same -> add gap in S2 
##		vertical transition -> j changes, i is the same -> add gap in S1  
##		diagonal transition -> i and j change -> no gaps

def findOptimalSequence():

	global datamatrix, M, N

	S1new = S2new = ""
	iCurr = iPrev = M
	jCurr = jPrev = N

	while (iCurr > 0 and jCurr > 0):

		# find the max cell location

		maxvalue = datamatrix[jCurr-1][iCurr-1]
		iMax = iCurr - 1
		jMax = jCurr - 1

		if (datamatrix[jCurr][iCurr-1] > maxvalue):
			maxvalue = datamatrix[jCurr][iCurr-1]
			iMax = iCurr - 1
			jMax = jCurr

		if (datamatrix[jCurr-1][iCurr] > maxvalue):
			maxvalue = datamatrix[jCurr-1][iCurr]
			iMax = iCurr
			jMax = jCurr - 1

		# place the correct character in the sequence

		if (iMax < iCurr and jMax < jCurr):
			S1new = S1[iMax] + S1new
			S2new = S2[jMax] + S2new

		elif (iMax < iCurr):
			S1new = S1[iMax] + S1new
			S2new = "-" + S2new

		elif (jMax < iCurr):
			S1new = "-" + S1new
			S2new = S2[jMax] + S2new	

		iCurr = iMax
		jCurr = jMax

	return S1new, S2new



if __name__ == '__main__':

	## sequences S1 and S2
	S1 = "CTAAGTACT"
	S2 = "CATTA"

	M = len(S1)		# num columns
	N = len(S2)		# num rows

	# the 2D matrix which stores the score for every sequence combo
	#		S1 ->
	# 	S2 	0	.	.
	# 	|	.	.	.
	# 	V 	.	.	.
	# must initialize datamatrix[0][0] = 0, otherwise everything else can 
	# be calculated
	datamatrix = [[0.0 for i in range(M+1)] for j in range(N+1)]

	## calculate the values of each cell in the datamatrix, top-left to bottom-right

	for i in range(0,N+1):			# O(M*N)

		for j in range(0,M+1):

			scoreCell(j,i)


	## find the sequence, bottom-right to top-left

	S1new, S2new = findOptimalSequence()


	## print results

	numDashes = max(len(S1), len(S2))

	print S1
	print S2
	print "-" * numDashes
	print S1new
	print S2new