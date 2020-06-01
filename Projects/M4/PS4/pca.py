import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, SparsePCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
from M2.PS2.kmeans import *

def readData(filename):

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
    
    return data, genes, patients


def pca(data, genes):

    dataMatrix = np.array(data)

    # scale data: center the average to 0
    # and scale so that standard deviation = 1
    scaledData = preprocessing.scale(dataMatrix)

    # do PCA
    pca = PCA(n_components=dataMatrix.shape[0])
    # fit uses SVD to create a model (U, SIGMA, V matrices) 
    # which only keeps n_components number of rows/cols
    pca.fit(scaledData)
    # transform actually applies the SVD transformation
    # to the data. it projects your data onto the n_component
    # number of PCs
    pcaData = pca.transform(scaledData)

    # extract results
    # calculate the % of variation accounted for by each PC
    perVariation = np.round(pca.explained_variance_ratio_*100, decimals=1)

    # bar graph
    labels = ["PC" + str(x) for x in range(1,len(perVariation)+1)]
    fig1 = plt.figure(1)
    plt.bar(x=range(1,len(perVariation)+1), height=perVariation, tick_label=labels)
    plt.ylabel("Percentage of explained variance")
    plt.xlabel("Principal component")
    plt.title("Scree plot")
    fig1.show()

    # scatter plot
    # fig2 = plt.figure(2)
    # plt.scatter(pcaData[:,0], pcaData[:,1])
    # plt.title("PCA graph")
    # plt.xlabel("PC1")
    # plt.ylabel("PC2")
    # fig2.show()

    pc2 = np.zeros(pcaData.shape[0], 2)
    pc2[:, 0] = pcaData[:, 0]
    pc2[:, 1] = pcaData[:, 1]
    kmeans(pc2, 2, "pca_plot")

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


def spca(data):

    return 0


if __name__ == "__main__":

    analysis_name = "tissue2"

    """creates directories for storing plots and intermediate files"""
    call(["rm", "-r", "./" + analysis_name + "_plots/"])
    call(["mkdir", "-p", "./" + analysis_name + "_plots/"])
    call(["mkdir", "-p", "./" + analysis_name + "_output/"])

    data, genes, patients = readData("M4/PS4/ExpData.txt")
    pcaData = pca(data, genes)

    

