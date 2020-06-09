import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, SparsePCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.abspath('../../M2/PS2'))
import kmeans

def readGenotypeData(filename):

    file = open(filename, 'r') 
    data = []
    patients = []

    # skip first line
    snps = file.readline()
    snps = snps.split("\t")
    snps.pop(0)

    while True:
        line = file.readline()
        if not line:
            break
        datastrings = line.split("\t")
        patients.append(datastrings.pop(0))
        datalist = [int(d) for d in datastrings]
        data.append(datalist)
    
    dataMatrix = np.array(data)
    
    return dataMatrix, snps, patients


def readExpressionData(filename):

    file = open(filename, 'r') 
    data = []
    patients = []

    # skip first line
    genes = file.readline()
    genes = genes.split("\t")
    genes.pop(0)

    while True:
        line = file.readline()
        if not line:
            break
        datastrings = line.split("\t")
        patients.append(datastrings.pop(0))
        datalist = [float(d) for d in datastrings]
        data.append(datalist)
    
    dataMatrix = np.array(data)
    
    return dataMatrix, genes, patients


def processResults(pca, pcaData, genes):

    # extract results
    # calculate the % of variation accounted for by each PC
    perVariation = np.round(pca.explained_variance_ratio_*100, decimals=1)

    # bar graph
    labels = ["PC" + str(x) for x in range(1,len(perVariation)+1)]
    # fig1 = plt.figure(1)
    # plt.bar(x=range(1,len(perVariation)+1), height=perVariation, tick_label=labels)
    # plt.ylabel("Percentage of explained variance")
    # plt.xlabel("Principal component")
    # plt.title("Scree plot")
    # fig1.show()

    # use k means to find clusters, k=2
    pc2 = np.zeros((pcaData.shape[0], 2))
    pc2[:, 0] = pcaData[:, 0]
    pc2[:, 1] = pcaData[:, 1]
    ptMem = kmeans.kmeans(pc2, 2, analysis_name, False)

    # scatter plot
    fig2 = plt.figure(2)
    pop1, pop2 = [], []
    for i in range(pcaData.shape[0]):
        if ptMem[i] == 0:
            pop1.append(pc2[i,:])
        elif ptMem[i] == 1:
            pop2.append(pc2[i,:])
    pop1arr = np.array(pop1)
    pop2arr = np.array(pop2)
    plt.scatter(pop1arr[:,0], pop1arr[:,1], c='b', marker="s", label='Pop 1')
    plt.scatter(pop2arr[:,0], pop2arr[:,1], c='r', marker="o", label='Pop 2')
    plt.title("PCA graph")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(loc='upper right');
    fig2.show()

    input()

    # find out which genes had the largest influence on 
    # separating the clusters along the PC1 axis
    # by looking at PC1 loading scores
    loadingScores = pd.Series(pca.components_[0], index=genes)
    loadingScores = loadingScores.abs().sort_values(ascending=False)
    top10 = loadingScores[0:10].index.values
    for i in range(10):
        try:
            print(top10[i], ": ", loadingScores[top10[i]])
        except ValueError:
            print(top10[i], " has multiple entries")


def preprocessData(data):

    # scale data: center the average to 0
    # and scale so that standard deviation = 1
    scaledData = preprocessing.scale(data)

    return scaledData


def pca(data, genes, analysis_name):

    preprocessedData = preprocessData(data)

    # do PCA
    pca = PCA(n_components=preprocessedData.shape[0])
    # fit uses SVD to create a model (U, SIGMA, V matrices) 
    # which only keeps n_components number of rows/cols
    pca.fit(preprocessedData)
    # transform actually applies the SVD transformation
    # to the data. it projects your data onto the n_component
    # number of PCs
    pcaData = pca.transform(preprocessedData)

    processResults(pca, pcaData, genes)
    

# too slow, too not use
# def spca(data, genes, analysis_name):

#     preprocessedData = preprocessData(data)

#     # do SPCA
#     spca = SparsePCA(n_components=10, alpha=1.0)
#     # fit uses SVD to create a model (U, SIGMA, V matrices) 
#     # which only keeps n_components number of rows/cols
#     spca.fit(preprocessedData)
#     # transform actually applies the SVD transformation
#     # to the data. it projects your data onto the n_component
#     # number of PCs
#     spcaData = spca.transform(preprocessedData)

#     processResults(spca, spcaData, genes)


if __name__ == "__main__":

    analysis_name = "pca"

    data, genes, patients = readExpressionData("ExpData.txt")
    pcaData = pca(data, genes, analysis_name)