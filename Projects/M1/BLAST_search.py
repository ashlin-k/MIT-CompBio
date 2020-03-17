import os
blosumdir = os.path.dirname(os.path.abspath(__file__)) + "\\blosum"
os.sys.path.insert(0,blosumdir) 
from blosum import *


def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):		# if string starts with > character, skip line
            continue
        seq.append(line.rstrip())		# else, append line to sequence and remove trailing characters, default is spaces
	
    return "".join(seq)					# converts array seq to string and returns


def getWords(seq, kmerlen):

	words = {}

	for i in range(0, len(seq) - kmerlen + 1):
		key = seq[i:i+kmerlen]
		words.setdefault(key, []).append(i)

	return words


def getNeighbourhood3(exactWord):

	global blosumMatrix

	threshold = 13

	nhWords = {}
	nhWords[exactWord] = 0
	for i in range(3):
		nhWords[exactWord] += blosumMatrix.lookup_score(exactWord[i], exactWord[i])

	allAminoAcids = "ARNDCQEGHILKMFPSTWYVBZX"
	testWord = "$$$"

	for i in range(0, len(exactWord)):

		testWord[i] = exactWord[i]

		for j in range(0, len(allAminoAcids)):

			if i == 0: index = 1
			else: index = 0
			testWord[index] = allAminoAcids[j]

			for k in range(0, len(allAminoAcids)): 

				if i == 2: index = 1
				else: index = 2
				testWord[index] = allAminoAcids[k]
				score = 0

				for l in range(0, len(exactWord)):
					score += blosumMatrix.lookup_score(testWord[l], exactWord[l])

				if score >= threshold:
					nhWords[testWord] = score

	return nhWords


def matchQueryWords(databaseWords, queryWords):

	iDatabaseMatches = []

	for qw in queryWords:
		if qw in databaseWords:
			iDatabaseMatches = iDatabaseMatches + databaseWords.get(qw)

	return iDatabaseMatches


def extendSequence(databaseString, queryString, kmerlen, databaseStart, queryStart):

	global blosumMatrix

	threshold = 0 # ? would think it would be lower than the other threshold

	iDatabase, iQuery = databaseStart, queryStart

	queryMatch = ""
	databaseMatch = ""

	# extend forward
	while iDatabase < len(databaseString) and iQuery < len(queryString):

		score += blosumMatrix.lookup_score(databaseString[iDatabase], queryString[iQuery])

		if score < threshold:
			return score, queryMatch, databaseMatch

		databaseMatch = databaseMatch + databaseMatch[iDatabase]

		if databaseString[iDatabase] == queryString[iQuery]:
			queryMatch = queryMatch + queryString[iQuery]
		else:
			queryMatch = queryMatch + "-"

		iDatabase++
		iQuery++


	# extend backward
	iDatabase, iQuery = databaseStart-1, queryStart-1
	while iDatabase >=0 and iQuery >= 0:

		score += blosumMatrix.lookup_score(databaseString[iDatabase], queryString[iQuery])

		if score < threshold:
			return score, queryMatch, databaseMatch

		databaseMatch = databaseMatch[iDatabase] + databaseMatch

		if databaseString[iDatabase] == queryString[iQuery]:
			queryMatch = queryString[iQuery] + queryMatch
		else:
			queryMatch = "-" + queryMatch

		iDatabase--
		iQuery--

	return score, queryMatch, databaseMatch

	

def main():

	global blosumMatrix
  
	# init blosum
	matrix_filename = blosumdir + "\\blosum62.txt"
	blosumMatrix = Matrix(matrix_filename)

	kmerlen = 3

	# generate database of words
	databaseFile = ""
	databaseString = readSeq(databaseFile)
	databaseWords = getWords(databaseString, kmerlen)

	# generate words from query
	queryFile = ""
	queryString = readSeq(queryFile)
	queryWords = getWords(queryString, kmerlen)

	# generate neighbourhood words
	for qw in queryWords:

		nhWords = getNeighbourhood3(qw)
		iQueryWords = queryWords.get(qw)
		iDatabaseMatches = matchQueryWords(databaseWords, nhWords)

		for i in range(0, len(iQueryWords)):

			for j in range(0, len(iDatabaseMatches)):
		
				queryStart = iQueryWords[i]
				databaseStart = iDatabaseMatches[j]
				score, queryMatch, databaseMatch = extendSequence(databaseString, queryString, kmerlen, databaseStart, queryStart)

				if score > threshold:
					# write to file or print


if __name__ == '__main__':
  main()
