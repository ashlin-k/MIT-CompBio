#!/usr/bin/env python

####
# 6.047/6.878 - Problem Set 1 - string hashing/dotplots
#
# INSTRUCTIONS FOR USE:
# call program as follows:
#  ./ps1-dotplot.py <FASTA 1> <FASTA 2> <PLOTFILE>
#     e.g. ./ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
#
# Make sure the ps1-dotplot.py is marked as executable:
#     chmod +x ps1-dotplot.py
# or in windows with:
#     python ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
# once you have put python in your path
#
#
# GNUPLOT
# Gnuplot is used to generate plots for this program.  It is a common plotting
# program installed on most unix systems.  To use gnuplot on athena do the
# following:
#
# athena% add gnu
#
# To test that it works do this:
#
# athena% gnuplot
# gnuplot> plot cos(x)
#
# you should then see a cosine plot appear.
#
# Note: plotting.py and util.py must be in the same directory as this script.
# These files contain the code for generating plots.  You should not
# worry about understanding any of the code contained within these files.
# Much of it is copied from another project and is unrelated.
#



import sys, random
import plotting


def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())

    return "".join(seq)


def quality(hits):
    """determines the quality of a list of hits"""

    # "good quality" means that the hits will be 
    # between two lines that intersect with (0,0)
    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    goodhits = []

    for hit in hits:
        upper = slope1 * hit[0] + offset1
        lower = slope2 * hit[0] + offset2

        if lower < hit[1] < upper:
            goodhits.append(hit)

    return goodhits


def makeDotplot(filename, hits):
    """generate a dotplot from a list of hits
       filename may end in the following file extensions:
         *.ps, *.png, *.jpg
    """

    if len(hits) == 0:
        return

    # zip splits the tuples of two in to a tuple of the
    # 0 index elements and a tuple of the 1 index elements 
    x, y = zip(* hits)

    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    hits2 = quality(hits)
    print "%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits)))

    # create plot
    p = plotting.Gnuplot()
    p.enableOutput(False)
    p.plot(x, y, xlab="sequence 2", ylab="sequence 1")
    p.plotfunc(lambda x: slope1 * x + offset1, 0, 1e6, 1e5)
    p.plotfunc(lambda x: slope2 * x + offset2, 0, 1e6, 1e5)

    # set plot labels
    p.set(xmin=0, xmax=1e6, ymin=0, ymax=1e6)
    p.set(main="dotplot (%d hits, %.5f%% hits on diagonal)" %
          (len(hits), 100 * len(hits2) / float(len(hits))))
    p.enableOutput(True)

    # output plot
    p.save(filename)

    return p


def invertSequence(seq):

    newSeq = ""

    for i in range(0, len(seq)):

        base = seq[i]
        if base == 'A':
            newSeq = 'T' + newSeq
        elif base == 'T':
            newSeq = 'A' + newSeq
        elif base == 'G':
            newSeq = 'C' + newSeq
        elif base == 'C':
            newSeq = 'G' + newSeq
    
    return newSeq


def findHits(seq1, seq2, kmerlen=30, matchspacing=1, inversion=False):

    # length of hash key

    # hash table for finding hits (dictionary)
    lookup = {}

    # store sequence hashes in hash table
    # if not in table, add empty list to the key and append index of loctiton
    # if already in table, just append index of location 
    print "hashing seq1..."
    for i in xrange(len(seq1) - kmerlen + 1):
        key = seq1[i:i+kmerlen:matchspacing]
        lookup.setdefault(key, []).append(i)

    # look up hashes in hash table
    print "hashing seq2..."
    hits = [] # an array of tuples (remember tuples are immutable)
    for i in xrange(len(seq2) - kmerlen + 1):

        if inversion:
            key = invertSequence(seq2[i:i+kmerlen:matchspacing])
        else:
            key = seq2[i:i+kmerlen:matchspacing]

        # store hits to hits list
        for hit in lookup.get(key, []):
            hits.append((i, hit))

    #
    # hits should be a list of tuples
    # [(index1_in_seq2, index1_in_seq1),
    #  (index2_in_seq2, index2_in_seq1),
    #  ...]
    #
    # for a one i in seq2, you could have multiple
    # hits in seq1
    #

    return hits


def main():

    # NOTE to WINDOWS users:
    #   If you do not want to use the command line, comment out the command line
    #   parsing code and hard-code the input filenames.
    #
    # For example, use the following:
    # file1 = "human-hoxa-region.fa"
    # file2 = "mouse-hoxa-region.fa"
    # plotfile = "dotplot.jpg"

    # parse command-line arguments
    if len(sys.argv) < 4:
        print "you must call program as:  "
        print "   python ps1-dotplot.py <FASTA 1> <FASTA 2> <PLOT FILE>"
        print "   PLOT FILE may be *.ps, *.png, *.jpg"
        sys.exit(1)
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    plotfile = sys.argv[3]



    # read sequences
    print "reading sequences"
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    # 2bi - 100, 1, False
    # 2bii - 60, 2, False
    # 2biii - 90, 3, False
    # 2biv - 120, 4, False
    # 2e - 180, 1, True -> 180 has only 1 hit
    kmerlen = 180 
    matchspacing = 1
    findInversion = True
    hits = findHits(seq1, seq2, kmerlen, matchspacing, findInversion)

    print "%d hits found" % len(hits)

    if (len(hits) > 0):
        print "making plot..."
        p = makeDotplot(plotfile, hits)


main()
