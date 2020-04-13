# Source for formulas:
# https://www.geeksforgeeks.org/ml-fuzzy-clustering/

import sys
import math
from subprocess import call
import numpy as np


def assignPoints(tbl, ctrs):
    """Assign each of the points in tbl to the cluster with
        center in ctrs"""

    global m

    # tbl is a Nx2 2D array, N is number of samples
    # ctrs is a kx2 2D array, where k = 3

    # I think ptsMembership is a Nxk 2D array, 
    # where ptsMembership[i] = [prob for k=0, prob for k=1, prob for k=2,...], where k = 0,1,2
    ptsMembership = []

    """SOME CODE GOES HERE"""
    
    K = len(ctrs)

    # this is the formula from lecture, but it doesn't work very well
    # p[k] = (1/(2*math.pi)) * math.exp(-dist/2)
    # so using formulas from above source

    for pt in tbl:

        # # get denom
        # Dsum = 0.0
        # for k in range(K):
        #     Dsum += euclideanDist(pt, ctrs[k])

        # if Dsum == 0 :
        #     print "Centers are all the same. Exiting function"
        #     exit(1)
            
        membership = [0.0 for i in range(K)]
        ptIsOnCtr = False

        # calculate num and membership
        for j in range(K):
            Dij = float(euclideanDist(pt, ctrs[j])) # num
            sumSq = 0.0
            for k in range(K):
                Dik = float(euclideanDist(pt, ctrs[k])) # denom
                if Dik != 0:
                    sumSq += (Dij**2) / (Dik**2)
                else:
                    membership = [0.0 for mem in range(K)]
                    membership[k] = 1.0
                    ptIsOnCtr = True
                    break
            
            if ptIsOnCtr:
                ptIsOnCtr = False
                break                
            
            membership[j] = 1.0 / math.pow(sumSq, 1/(m-1))
        
        sumMem = sum(membership)
        ptsMembership.append(membership)
      
    return ptsMembership


def recalculateCtrs(tbl, ctrs, ptsMembership):
    """Update the centroids based on the points assigned to them"""

    global m

    K = len(ctrs)
    D = len(ctrs[0])
    newCtrs = [[0.0 for d in range(D)] for k in range(K)]
    pSum = [0.0 for k in range(K)]

    """SOME CODE GOES HERE"""
    
    # sum
    for i in range(len(ptsMembership)):
        for k in range(K):
            pSum[k] += ptsMembership[i][k]**m
            for d in range(D):
                newCtrs[k][d] += tbl[i][d] * ptsMembership[i][k]**m
    
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


def getRandomCtr(limits):

    newCtr = []

    # assuming n features, limits is organized as:
    # [
    #   [x1min, x1max],
    #   [x2min, x2max],
    #   ... ,
    #   [xnmin, xnmax]
    # ]

    for i in range(len(limits)):
        xMin = limits[i][0]
        xMax = limits[i][1]
        newCtr.append(np.random.uniform(xMin, xMax))

    return newCtr


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

    # m is a FKM param normally set to 2 when the distribution is not known
    # 1 < m < inf
    global m    
    m = 2       

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
    K = 3
    minX1, minX2, maxX1, maxX2 = 0, 0, 0, 0
    for i in range(len(dataTable)):
        if dataTable[i][0] < minX1:
            minX1 = dataTable[i][0]
        elif dataTable[i][0] > maxX1:
            maxX1 = dataTable[i][0]
        if dataTable[i][1] < minX2:
            minX2 = dataTable[i][1]
        elif dataTable[i][1] > maxX2:
            maxX2 = dataTable[i][1]
    limits = [[minX1, maxX1], [minX2, maxX2]]
    newCtrs = []
    for k in range(K):
        newCtrs.append(getRandomCtr(limits))
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


