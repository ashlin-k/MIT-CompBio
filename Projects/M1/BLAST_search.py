import os
blosumdir = os.path.dirname(os.path.abspath(__file__)) + "\\blosum"
os.sys.path.insert(0,blosumdir) 
from blosum import *


class Match:

	def __init__(self, iDB, iQ, nDB, nQ, dbString="", qString="", score=0):
		self._iDatabase = iDB
		self._iQuery = iQ
		self._Ndatabase = nDB
		self._Nquery = nQ
		self._databaseString = dbString
		self._queryString = qString
		self._score = score


	def setDatabaseSequence(self, seq):
		self._databaseString = seq

	def setQuerySequency(self, seq):
		self._queryString = seq

	def setScore(self, score):
		self._score = score

	def getDatabaseSequence(self):
		return self._databaseString

	def getDatabaseIndex(self):
		return self._iDatabase

	def getQuerySequence(self):
		return self._queryString

	def getQueryIndex(self):
		return self._iQuery

	def getScore(self):
		return self._score


class Hits:

	def __init__(self):
		self._matches = {}	# keys will be the DB index

	def hasMatch(self, match):
		return match.getDatabaseIndex() in self._matches

	def addMatch(self, match, kmerlen):
		
		# if match is already in the database, don't add
		if self.hasMatch(match):
			return False

		else:
			iNew = match.getDatabaseIndex()
			for i in range(iNew - kmerlen, iNew + kmerlen + 1):

				# if adjacent duplicate is found, don't add
				if i in self._matches:
					return False

			# no adjacent duplicate is found, therefore add
			self._matches[iNew] = match
			return True

	def getHitDictionary(self):
		return self._matches
			


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

	thresholdT = 13

	nhWords = {}
	exactWordString = "".join(exactWord)
	nhWords[exactWordString] = 0
	for i in range(3):
		nhWords[exactWordString] += int(blosumMatrix.lookup_score(exactWord[i], exactWord[i]))

	allAminoAcids = "ARNDCQEGHILKMFPSTWYVBZX"
	testWordString = "$$$"
	testWord = list(testWordString)

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
					score += int(blosumMatrix.lookup_score(testWord[l], exactWord[l]))

				if score >= thresholdT:
					nhWords[testWordString] = score

	return nhWords


def matchQueryWords(databaseWords, queryWords):

	iDatabaseMatches = []

	for qw in queryWords:
		if qw in databaseWords:
			iDatabaseMatches = iDatabaseMatches + databaseWords.get(qw)

	return iDatabaseMatches


def extendSequence(databaseString, queryString, kmerlen, databaseStart, queryStart):

	global blosumMatrix, thresholdX

	iDatabase, iQuery = databaseStart, queryStart

	queryMatch = ""
	databaseMatch = ""
	score = 0

	# extend forward
	while iDatabase < len(databaseString) and iQuery < len(queryString):

		score += int(blosumMatrix.lookup_score(databaseString[iDatabase], queryString[iQuery]))

		if score < thresholdX:
			match = Match(databaseStart, queryStart, len(databaseMatch), len(queryMatch), databaseMatch, queryMatch, score)
			return match

		databaseMatch = databaseMatch + databaseString[iDatabase]

		if databaseString[iDatabase] == queryString[iQuery]:
			queryMatch = queryMatch + queryString[iQuery]
		else:
			queryMatch = queryMatch + "-"

		iDatabase += 1
		iQuery += 1


	# extend backward
	iDatabase, iQuery = databaseStart-1, queryStart-1
	while iDatabase >=0 and iQuery >= 0:

		score += int(blosumMatrix.lookup_score(databaseString[iDatabase], queryString[iQuery]))

		if score < thresholdX:
			match = Match(iDatabase, iQuery, len(databaseMatch), len(queryMatch), databaseMatch, queryMatch, score)
			return match

		databaseMatch = databaseString[iDatabase] + databaseMatch

		if databaseString[iDatabase] == queryString[iQuery]:
			queryMatch = queryString[iQuery] + queryMatch
		else:
			queryMatch = "-" + queryMatch

		iDatabase -= 1
		iQuery -= 1


	match = Match(iDatabase+1, iQuery+1, len(databaseMatch), len(queryMatch), databaseMatch, queryMatch, score)
	return match

	

def main():

	global blosumMatrix, thresholdX
  
	# init blosum
	matrix_filename = blosumdir + "\\blosum62.txt"
	blosumMatrix = Matrix(matrix_filename)

	kmerlen = 3

	# generate database of words
	databaseFile = "database.fa"
	databaseString = readSeq(databaseFile)
	# databaseString = "TPQRRABPVDWWRABPT"
	databaseWords = getWords(databaseString, kmerlen)

	# generate words from query
	queryFile = "query.fa"
	queryString = readSeq(queryFile)
	# queryString = "VDW"
	queryWords = getWords(queryString, kmerlen)

	thresholdX = 0

	# generate neighbourhood words
	# find seeds and extend
	hits = Hits()
	for qw in queryWords:

		nhWords = getNeighbourhood3(list(qw))
		iQuery = queryWords.get(qw)		# a list of all locations where the word is found in the query
										# assumes you may find the same k-mer more than once in a given query, 
										# but probably it will be 1

		# find seeds
		iDatabaseMatches = matchQueryWords(databaseWords, nhWords)

		# extend
		for i in range(0, len(iQuery)):	

			for j in range(0, len(iDatabaseMatches)):
		
				queryStart = iQuery[i]
				databaseStart = iDatabaseMatches[j]
				match = extendSequence(databaseString, queryString, kmerlen, databaseStart, queryStart)

				if match.getScore() > thresholdX:
					hits.addMatch(match, kmerlen)


	# print all matches
	hitDictionary = hits.getHitDictionary()
	for hit in hitDictionary:
		print "iDatabase = ", hit, ", iQuery = ", hitDictionary[hit].getQueryIndex(), ", Score = ", hitDictionary[hit].getScore()
		print "DB: ", hitDictionary[hit].getDatabaseSequence()
		print "QY: ", hitDictionary[hit].getQuerySequence()
		print "-" * 15
					


if __name__ == '__main__':
	main()
