## PS4 problem 1: suffix bwt algorithm / Burrows Wheeler Transform (BWT)
## The aim of BWT is to map millions of short reads to the genome in O(m) time
## where m is the length of the short read

from SuffixTree import SuffixTree

if __name__ == "__main__":

    genome = "MANOLISKELLIS"
    bwt = SuffixTree(genome)
    sp, ep, bwtsp, bwtep = bwt.findMatch("OLAS") 

    if sp > 0 and ep > 0:

        print("Table ptr locations:")
        print ("sp = ", sp, ", ep = ", ep)
        print("-----")   

        locations = bwt.getGenomeLocations(bwtsp, bwtep)

        print("Genome locations:")
        print ("bwtsp = ", bwtsp, ", bwtep = ", bwtep)
        print(genome)
        i, j = 0, 0
        while j < len(genome) and i < len(locations):
            if locations[i] == j:
                print("^", end = '')
                i += 1
            else:
                print(" ", end = '')
            j += 1
        print("")
    

