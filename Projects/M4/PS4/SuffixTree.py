class SuffixTree:

    def __init__(self, genome=""):

        self.bwt = []
        self.bwtOrdered = []
        self.S = []
        self.isInit = False
        self.markerChar = "$"

        if genome != "":
            self.initBwt(genome)
    

    def initBwt(self, genome):

        if self.isInit:
            print("This SuffixTree has already been init")
            return self.bwt

        print("Genome = ", genome)

        augGenome = genome + self.markerChar
        M = len(augGenome)
        table = []
        
        # create tree
        for i in range(M):
            table.append(augGenome)
            augGenome = augGenome[M-1] + augGenome[0:-1]
        
        # sort rows alphabetically
        self.bwtSort(table)
        print("-----")
        for i in range(len(table)):
            loc = len(table[i]) - table[i].find(self.markerChar) - 2
            self.S.append(loc)
            print(i, "\t", table[i])
        print("-----")
        
        # take last row and store it as an array
        self.bwt = [row[-1:] for row in table]
        print(self.bwt)
        print(self.S)
        print("-----")

        # create an ordered array for fn C(a)
        self.bwtOrdered = self.bwt.copy()
        self.bwtSort(self.bwtOrdered)

        self.isInit = True

        return self.bwt
    

    def bwtSort(self, tree):

        tree.sort()
        originalLength = len(tree)
        tree.insert(0,"")

        for i in range(1, originalLength+1):

            if tree[i][0] == self.markerChar:   # last
                tree[0] = tree.pop(i)
            
            if len(tree) == originalLength:
                break

        return tree


    # Let C(a): number of characters lexicographically smaller than a
    # in the string.
    # For example, given the word BANANA$, C($) = 0, C(A) = 1, C(B) = 4,
    # and C(N) = 5.
    def C(self, a, n=0):

        index = self.bwtOrdered.index(a)

        if index == len(self.bwtOrdered) - 1:
            return index

        if n != 0:
            a2 = a
            while index < len(self.bwtOrdered)-1:        
                index += 1        
                a2 = self.bwtOrdered[index]
                if a2 != a:
                    return index - 1               
        
        return index

    # Returns of the number of occurences of char c
    # before index k in the bwt array
    def occurrence(self, c, k):

        count = 0
        for i in range(k):
            if self.bwt[i] == c:
                count += 1
        
        return count


    def findExactMatch(self, motif):

        if not self.isInit:
            return -1,-1,-1,-1

        if len(motif) <= 1:
            print("Motif must be at least two chars in length")
            return -1,-1,-1,-1

        N = len(self.bwt)
        M = len(motif)
        
        c = motif[M-1]
        if c not in self.bwt:
            print("Motif contains invalid charachter")
            return -1,-1,-1,-1
        sp = self.C(c)      # index of first appearance of last char
        ep = self.C(c, 1)   # index of first appearance of next lexographical char
        bwtsp = sp
        bwtep = ep

        i = M - 2;
        # print(self.bwt)
        # print(self.bwtOrdered)
        # print("c=",c ,", sp=",sp, ", ep=", ep)

        while sp <= ep and i >= 0:
            c = motif[i]
            if c not in self.bwt:
                print("Motif contains invalid charachter")
                return -1,-1,-1,-1
            # print("c=",c , ", C(c)=", self.C(c), \
            #     ", occ(c,sp)=",self.occurrence(c, sp), \
            #     ", occ(c,ep)=",self.occurrence(c, ep+1), end="")     
            bwtsp = sp
            bwtep = ep       
            sp = self.C(c) + self.occurrence(c, sp)
            ep = self.C(c) + self.occurrence(c, ep+1) - 1
            # print(", sp=",sp, ", ep=", ep)
            i -= 1
        
        if sp == ep:
            if c == self.bwt[bwtsp]:
                bwtep = bwtsp
            else:
                bwtsp = bwtep

        return sp, ep, bwtsp, bwtep

    
    def ibwt(self):

        if not self.isInit:
            return ""
        
        N = len(self.bwt)
        table = ["" for i in range(N)]

        for i in range(N):

            # add last col array
            for j in range(N):
                table[j] = self.bwt[j] + table[j]
            
            # sort
            self.bwtSort(table)
        
        # will be the first entry since table is sorted
        # and marker char is lexographically first
        genome = table[0].replace(self.markerChar, "")

        return genome
    

    def getGenomeLocations(self, bwtsp, bwtep):

        locations = []
        i = bwtsp
        while i <= bwtep:
            locations.append(self.S[i])
            i += 1

        locations.sort()
        return locations
    

    