# Support vector machine
# This program uses data from PS2/tissue1_data.txt training
# The kernel function here is the RBF

from LinearSVM import LinearSVM
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets.samples_generator import make_blobs
from sklearn.model_selection import train_test_split
from sklearn.svm import LinearSVC
from sklearn.metrics import confusion_matrix


def f(x, w, b, c=0):
    return (-w[0] * x - b + c) / w[1]


def getDataFromFile(filename):

    """Reads in the point data from the given tissue file"""
    dataTable = []
    f = open("./" + analysis_name + "_data.txt")
    for dataLine in f:
        dataTable.append([float(str) for str in dataLine.rstrip().split("\t")])
    f.close()

    return dataTable


def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    # if len(sys.argv) < 1:
    #     print "you must call program as: python ./kmeans.py <datafile>"
    #     sys.exit(1)
    # analysis_name = sys.argv[1]
    # analysis_name = "tissue2"

    # generate fake data points
    X, y = make_blobs(n_samples=250, centers=2,
                  random_state=0, cluster_std=0.60)
    y[y == 0] = -1  # if y == 0, set y = -1 instead
    tmp = np.ones(len(X))
    y = tmp * y
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # create and train the svm
    svm = LinearSVM()
    svm.train(X_train, y_train)

    # plot the data
    plt.scatter(X_train[:, 0], X_train[:, 1], c=y_train, cmap='winter')
    # w.x + b = 0 line
    a0 = -4; a1 = f(a0, svm.w, svm.b)
    b0 = 4; b1 = f(b0, svm.w, svm.b)
    plt.plot([a0,b0], [a1,b1], 'k')
    # w.x + b = 1 line
    a0 = -4; a1 = f(a0, svm.w, svm.b, 1)
    b0 = 4; b1 = f(b0, svm.w, svm.b, 1)
    plt.plot([a0,b0], [a1,b1], 'k--')
    # w.x + b = -1 line
    a0 = -4; a1 = f(a0, svm.w, svm.b, -1)
    b0 = 4; b1 = f(b0, svm.w, svm.b, -1)
    plt.plot([a0,b0], [a1,b1], 'k--')
    plt.show()

    # do a test sample
    y_pred = svm.predict(X_test)

    # get the confusion matrix, which tells us the number of
    # predicted (rows) vs actual (cols) instances of classes
    conf = confusion_matrix(y_test, y_pred)
    print(conf)
    


if __name__ == "__main__":
    main()