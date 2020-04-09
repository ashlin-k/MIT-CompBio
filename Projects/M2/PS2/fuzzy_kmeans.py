import sys
import math
from subprocess import call
import numpy as np


def assignPoints(tbl, ctrs):
    """Assign each of the points in tbl to the cluster with
        center in ctrs"""

    # tbl is a Nx2 2D array, N is number of samples
    # ctrs is a kx2 2D array, where k = 3

    # I think ptsAsgn is a Nxk 2D array, 
    # where ptsAsgn[i] = [prob for k=0, prob for k=1, prob for k=2,...], where k = 0,1,2
    ptsAsgn = []

    """SOME CODE GOES HERE"""
    
    K = len(ctrs)

    # this is the formula from lecture, but it doesn't work very well
    # p[k] = (1/(2*math.pi)) * math.exp(-dist/2)
    # so using formula from wikipedia https://en.wikipedia.org/wiki/Fuzzy_clustering#Algorithm

    m = 2       # m is a param normally set to 2 when the distribution is not known

    for pt in tbl:

        # get denom
        denom = 0.0
        for k in range(K):
            denom += euclideanDist(pt, ctrs[k])

        weight = [0 for i in range(K)]

        # calculate num and weight
        for k in range(K):
            num = float(euclideanDist(pt, ctrs[k]))
            denom2 = math.pow((num/denom), 2*m-1)
            if denom2 != 0:
                weight[k] = 1 / math.pow((num/denom), 2*m-1)
            else:       # if pt == ctrs[k], set that probability to 1 and rest to zero
                weight = [0 for i in range(K)]
                weight[k] = 1

        ptsAsgn.append(weight)

    return ptsAsgn


def recalculateCtrs(tbl, ctrs, ptsAsgn):
    """Update the centroids based on the points assigned to them"""

    K = len(ctrs)
    D = len(ctrs[0])
    newCtrs = [[0.0 for d in range(D)] for k in range(K)]
    pSum = [0.0 for k in range(K)]

    """SOME CODE GOES HERE"""
    
    # sum
    for i in range(len(ptsAsgn)):
        for k in range(K):
            pSum[k] += ptsAsgn[i][k]
            for d in range(D):
                newCtrs[k][d] += tbl[i][d] * ptsAsgn[i][k]
    
    # divide by total
    for k in range(K):
        for d in range(D):
            newCtrs[k][d] /= pSum[k]

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

    p = open("./fuzzy_" + anLabel + "_output/dummy_table.txt", "w")

    for i in xrange(len(tbl)):
        for j in xrange(len(tbl[i])):
            p.write(`tbl[i][j]`)
            p.write("\t")
        k = ptMemb[i].index(max(ptMemb[i]))
        p.write(`k`)
        p.write("\n")

    for i in xrange(len(cntrs)):
        for j in xrange(len(cntrs[i])):
            p.write(`cntrs[i][j]`);
            p.write("\t")
        p.write("Clust" + `i`)
        p.write("\n")

    p.close()

    plotCMD = "R CMD BATCH '--args ./fuzzy_" + anLabel + "_output/dummy_table.txt ./fuzzy_" + anLabel \
        + "_plots/cluster_step%d.png" % stepCnt + "' ./kmeans_plot.R;"
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
    call(["rm", "-r", "./fuzzy_" + analysis_name + "_plots/"])
    call(["mkdir", "-p", "./fuzzy_" + analysis_name + "_plots/"])
    call(["mkdir", "-p", "./fuzzy_" + analysis_name + "_output/"])

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
        if stopDist < 1:
            stopCrit = True

        stepCount = stepCount + 1

if __name__ == "__main__":
    main()


