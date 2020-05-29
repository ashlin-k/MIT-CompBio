import numpy as np
from sklearn.decomposition import PCA

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



if __name__ == "__main__":

    data, genes, patients = readData("M4/PS4/ExpData.txt")
    
    pca = PCA(n_components=data.shape[0])
    pca.fit(data)