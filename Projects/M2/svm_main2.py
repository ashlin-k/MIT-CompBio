# Support vector machine using the python library sklearn
# Tutorial linear: https://www.datacamp.com/community/tutorials/svm-classification-scikit-learn-python
# Tutorial RBF: https://chrisalbon.com/machine_learning/support_vector_machines/svc_parameters_using_rbf_kernel/,
#   https://scikit-learn.org/stable/auto_examples/svm/plot_rbf_parameters.html

# SVM tuning params
# C: Trades off correct classification of training examples against 
#   maximization of the decision fn's margin. High C is more correct for
#   more complex decision function, low C is less correct but simpler decision fn
# Gamma: defines how far the influence of a single training example reaches; 
#   can be seen as the inverse of the radius of influence of samples selected 
#   by the model as support vectors

# One thing we didn't do in rbfSVM, but should is normalize the training data
# using sklearn.preprocessing.StandardScaler:
#   scalar = StandardScaler()
#   X = scalar.fit_transform(X)
# It scales data by subtracting the mean and dividing by the standard deviation. 
# Therefore it focuses on the distribution of points rather than on the range. 
# You can then scale the test data back afterwards.
# https://stats.stackexchange.com/questions/65094/why-scaling-is-important-for-the-linear-svm-classification


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap

from sklearn import svm, datasets, metrics
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC


def versiontuple(v):
    return tuple(map(int, (v.split("."))))


def plot_decision_regions(X_train, X_test, y_train, y_test, classifier, test_idx=None, resolution=0.02):

    # setup marker generator and color map
    markers = ('s', 'x', 'o', '^', 'v')
    colors = ('red', 'blue', 'lightgreen', 'gray', 'cyan')
    cmap = ListedColormap(colors[:len(np.unique(y_train))])

    # plot the decision surface
    x1_min, x1_max = X_train[:, 0].min() - 1, X_train[:, 0].max() + 1
    x2_min, x2_max = X_train[:, 1].min() - 1, X_train[:, 1].max() + 1
    xx1, xx2 = np.meshgrid(np.arange(x1_min, x1_max, resolution),
                           np.arange(x2_min, x2_max, resolution))
    Z = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
    Z = Z.reshape(xx1.shape)
    plt.contourf(xx1, xx2, Z, alpha=0.4, cmap=cmap)
    plt.xlim(xx1.min(), xx1.max())
    plt.ylim(xx2.min(), xx2.max())

    for idx, cl in enumerate(np.unique(y_train)):
        labelTrain = "y = " + str(cl) + ", train"
        labelTest = "y = " + str(cl) + ", test"
        plt.scatter(x=X_train[y_train == cl, 0], y=X_train[y_train == cl, 1], \
            alpha=0.8, c=cmap(idx), marker=markers[2*idx], label=labelTrain)
        plt.scatter(x=X_test[y_test == cl, 0], y=X_test[y_test == cl, 1], \
            alpha=0.8, c=cmap(idx), marker=markers[2*idx+1], label=labelTest)

    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.show()


def evaluateSVM(y_test, y_pred):

    # check accuracy
    # precision is what % of positive tuples are labeled as positive
    acc = metrics.accuracy_score(y_test, y_pred)
    prec = metrics.precision_score(y_test, y_pred)
    print("Accuracy = ", acc)
    print("Precision = ", prec)


def linearSVM():

    # load dataset
    cancerData = datasets.load_breast_cancer()

    # split the data into training and testing data
    # 70% training and 30% test
    X_train, X_test, y_train, y_test = \
        train_test_split(cancerData.data, cancerData.target, test_size=0.3,random_state=109) 

    # svm classifier
    clf = svm.SVC(kernel='linear')

    # train the model
    clf.fit(X_train, y_train)

    # predict the response from test datastet
    y_pred = clf.predict(X_test)

    # evaluate accuracy
    evaluateSVM(y_test, y_pred)


def rbfSVM():

    # generate fake data

    np.random.seed(0)
    X_xor = np.random.randn(200, 2)
    y_xor = np.logical_xor(X_xor[:, 0] > 0,
                        X_xor[:, 1] > 0)
    y_xor = np.where(y_xor, 1, -1)

    # split the data into training and testing data
    # 70% training and 30% test
    X_train, X_test, y_train, y_test = \
        train_test_split(X_xor, y_xor, test_size=0.3,random_state=109) 

    # find best C and gamma parms

    C_range = np.logspace(-2, 10, 13)
    gamma_range = np.logspace(-9, 3, 13)
    param_grid = dict(gamma=gamma_range, C=C_range)
    cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv)
    grid.fit(X_train, y_train)

    print("The best parameters are %s with a score of %0.2f"
        % (grid.best_params_, grid.best_score_))

    # create classifier

    C = grid.best_params_['C']
    gamma = grid.best_params_['gamma']
    svm = SVC(kernel='rbf', random_state=0, gamma=gamma, C=C)
    svm.fit(X_train, y_train)
    y_pred = svm.predict(X_test)
    plot_decision_regions(X_train, X_test, y_train, y_pred, classifier=svm)

    # evaluate accuracy
    evaluateSVM(y_test, y_pred)


    



if __name__ == "__main__":
    rbfSVM()