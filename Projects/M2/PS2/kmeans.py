import sys
import math
from subprocess import call
import numpy as np


def assignPoints(tbl, ctrs):
    """Assign each of the points in tbl to the cluster with
        center in ctrs"""

    # tbl is a Nx2 2D array, N is number of samples
    # ctrs is a kx2 2D array, where k = 3

    # I think ptsAsgn is a 1D array, where ptsAsign[i] = k, where k = 0,1,2
    ptsAsgn = []

    """SOME CODE GOES HERE"""
    
    K = len(ctrs)

    for pt in tbl:
        dist = [0.0 for d in range(K)]
        for k in range(K):
            dist[k] = euclideanDist(pt, ctrs[k])
        ptsAsgn.append(dist.index(min(dist)))

    return ptsAsgn


def recalculateCtrs(tbl, ctrs, ptsAsgn):
    """Update the centroids based on the points assigned to them"""

    K = len(ctrs)
    D = len(ctrs[0])
    newCtrs = [[0.0 for j in range(D)] for k in range(K)]

    """SOME CODE GOES HERE"""

    counts = [0 for k in range(K)]
    
    # sum
    for i in range(len(ptsAsgn)):
        counts[ptsAsgn[i]] += 1
        for d in range(D):
            newCtrs[ptsAsgn[i]][d] += tbl[i][d]
    
    # divide by total
    for k in range(K):
        for d in range(D):
            newCtrs[k][d] /= counts[k]

    return newCtrs


def euclideanDist(x, y):

    if len(x) != len(y):
        print "x and y must be the same length"
        sys.exit(1)
    else:
        dist_val = 0
        for i in xrange(len(x)):
            dist_val = dist_val + math.pow((x[i] - y[i]),2)

    return math.sqrt(dist_val)


def plotClusters(tbl, ptMemb, cntrs, stepCnt, anLabel):
    """generate a scatterplot of the current
       k-means cluster assignments

       filename may end in the following file extensions:
         *.ps, *.png, *.jpg
    """

    p = open("./" + anLabel + "_output/dummy_table.txt", "w")

    for i in xrange(len(tbl)):
        for j in xrange(len(tbl[i])):
            p.write(`tbl[i][j]`)
            p.write("\t")
        p.write(`ptMemb[i]`)
        p.write("\n")

    for i in xrange(len(cntrs)):
        for j in xrange(len(cntrs[i])):
            p.write(`cntrs[i][j]`);
            p.write("\t")
        p.write("Clust" + `i`)
        p.write("\n")

    p.close()

    plotCMD = "R CMD BATCH '--args ./" + anLabel + "_output/dummy_table.txt ./" + anLabel + "_plots/cluster_step%d.png" % stepCnt + "' ./kmeans_plot.R;"
    call(plotCMD, shell=True)



###############################################################################
# MAIN
###############################################################################

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    # if len(sys.argv) < 1:
    #     print "you must call program as: python ./kmeans.py <datafile>"
    #     sys.exit(1)
    # analysis_name = sys.argv[1]
    analysis_name = "tissue2"

    """creates directories for storing plots and intermediate files"""
    call(["rm", "-r", "./" + analysis_name + "_plots/"])
    call(["mkdir", "-p", "./" + analysis_name + "_plots/"])
    call(["mkdir", "-p", "./" + analysis_name + "_output/"])

    """Reads in the point data from the given tissue file"""
    dataTable = []
    f = open("./" + analysis_name + "_data.txt")
    for dataLine in f:
        dataTable.append([float(str) for str in dataLine.rstrip().split("\t")])
    f.close()

    """initializes centroids, stop criterion and step counting for clustering"""
    newCtrs = []
    for i in range(3):
        rdm = np.random.random_integers(0,len(dataTable)-1)
        newCtrs.append(dataTable[rdm])
    # newCtrs = [[5,0], [5,40], [5,80]]
    ptMemb = assignPoints(dataTable, newCtrs)
    stopCrit = False
    stepCount = 1

    """performs k-means clustering, plotting the clusters at each step"""
    while stopCrit == False:
        plotClusters(dataTable, ptMemb, newCtrs, stepCount, analysis_name)

        """SOME CODE GOES HERE"""
        
        oldCtrs = newCtrs
        newCtrs = recalculateCtrs(dataTable, newCtrs, ptMemb)
        ptMemb = assignPoints(dataTable, newCtrs)


        """stop criterion - when centroids' total movement after a step is below
            the threshold, stop the algorithm"""
        stopDist = 0
        for i in xrange(len(newCtrs)):
            stopDist = stopDist + euclideanDist(oldCtrs[i], newCtrs[i])
        if stopDist < 5:
            stopCrit = True

        stepCount = stepCount + 1

if __name__ == "__main__":
    main()


