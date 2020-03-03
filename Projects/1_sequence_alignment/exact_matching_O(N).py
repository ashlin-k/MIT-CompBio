
#	This program will find an exact number of matches of an M-length
#	pattern, P, in a N-length string, T, in linear time, O(N) using the
# 	Z algorithm, Boyer-Moore and Knuth-Morris-Pratt (KMP)

#	Based of of Lecture 3 in the MIT CompBio lecture series
#	https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-047-computational-biology-fall-2015/lectures_slides/MIT6_047F15_Lecture03.pdf

# 	ZString:
#	A string that can be preprocessed to store a Z array.
# 	Z[i] = length of longest substirng starting from Z[i],
# 	which is also a prefix of the string
# 	(meaning the first char of the substring must be str[0])

class ZString:

	def __init__(self, str):
		self.string = str 						# the string
		self.N = len(self.string)				# size of the string
		self.Z = [0 for i in range(self.N)]		# Z values, redundancy structure, init to 0
		self.ZisInit = False					

	def size(self):
		return self.N

	def getString(self):
		return self.string

	def getZ(self):
		return self.Z

	# code borrowed from GeeksForGeeks
	# https://www.geeksforgeeks.org/z-algorithm-linear-time-pattern-searching-algorithm/
	def preprocessZ(self):

		l, r, k = 0, 0, 0

		for i in range(1, self.N): 

			if i > r: 
				l, r = i, i 
				while r < self.N and self.string[r - l] == self.string[r]: 
					r += 1
				self.Z[i] = r - l
				r -= 1

	    	else: 

	    		k = i - l 

	    		if self.Z[k] < r - i + 1:
	    			self.Z[i] = self.Z[k] 

	    		else: 
	    			l = i 
	    			while r < self.N and self.string[r - l] == self.string[r]: 
	    				r += 1
	    			self.Z[i] = r - l 
	    			r -= 1

		self.ZisInit = True

		return self.Z




if __name__ == '__main__':

	# the pattern, P, and the text, T
	P = "ACA"					# the pattern string, size M
	T = "TAACAGACACAGG"			# the text string, size N

	# preprocess concatenated P-T string, O((M+N)^2)
	concatenatedPT = P + "$" + T
	PT = ZString(concatenatedPT)
	z = PT.preprocessZ()

	# find P in T, O(M+N)
	M = len(P)
	indicatorString = ""
	arrowString = "^" * M
	startArrow = M-3
	for i in range(M+1, PT.size()):

		if z[i] == M:
			indicatorString += arrowString
			startArrow = i

		else:
			if i > startArrow + 2:
				indicatorString += " "


	print "P = ", P
	print "T = ", T
	print " " * len("T = "), indicatorString
			

	

