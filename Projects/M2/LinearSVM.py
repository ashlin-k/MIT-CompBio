# Support Vector Machine class

# reference: https://towardsdatascience.com/support-vector-machine-python-example-d67d9b63f1c8

# We will frame this as a max/min problem. From the lectures we know that we must
#   max 2/|w| -> min|w| -> min 0.5*|w|^2 (take the integral)

# Next we can use the Lagrangian multiplier to solve the min/max problem.
# https://en.wikipedia.org/wiki/Lagrange_multiplier
# Recall that for a cost function, f(x), with constraint g(x), the 
# Lagrangian is:
#   L = f(x) - a*g(x),  where a is a multiplier
# We can then calculate the partial derivatives of L and set them to 0
# to find a local max/min

# For this problem, 
#   f(w) = 0.5|w|^2
#   g(w, b) = sum(i=1 to n) [(yi*(w.x + b) - 1)]
#   a = ai, i = 1 to n
# Therefore
#   L = 0.5|w|^2 - sum(i=1 to n) [ai * (yi*(w.x + b) - 1)]
# Given the identities
#   w = sum(i=1 to n) [ai*yi*xi]    (1)
#   sum(i=1 to n) [ai*yi] = 0       (2)
# We can simplify L (in matrix form) to
#   L = 1.' * a - 0.5 * a.' * X.' * X * a
#   for constraints y.' * a = 0 and a >= 0
# where 
#   A.' is the transpose of matrix A
#   1 = [1 1 ... 1].'               (Nx1)
#   a = [a1 a2 ... an].'            (Nx1)
#   X = [x1*y1 x2*y2 ... xn*yn]     (fxN)
#   xi is a (fx1) vector
#   yi is a scalar = 1 or -1

# We can use the CVXOPT library to solve for max a
# We need to provide it with the parameters for this type of eqn:
#   min(x) [q.' * x + 0.5 * x.' * P * x]
#   for constraints A*x = b and G*x <= h
# Our equivalent parameters will maximize a instead of minimize x,
# and equivalent P and q are:
#   P => X.' * X
#   q => -1 = [-1]  (Nx1)
#   A => y.'        (1xN)
#   b => 0          (1x1)
#   G => diag[-1]   (NxN)
#   h => 0 = [0]    (Nx1)

import numpy as np
import cvxopt

class LinearSVM:

    def train(self, X, y):

        # X is a Nxf matrix (f = numFeatures, N = numSamples)
        # y is a Nx1 matrix

        numSamples, numFeatures = X.shape   # shape returns rows, cols

        # P => X.' * X
        # where X = [x1*y1 x2*y2 ... xn*yn].' => Xy  
        XyT = np.zeros((numSamples, numFeatures))
        for i in range(numSamples):
            XyT[i] = X[i] * y[i]
        Xy = XyT.transpose()
        P = cvxopt.matrix(XyT.dot(Xy))

        # q => -1 = [-1]  (Nx1)
        q = cvxopt.matrix(np.ones(numSamples) * -1)

        # A => y.'      (1xN)
        A = cvxopt.matrix(y, (1,numSamples))

        # b => 0 = [0]    (1x1)
        b = cvxopt.matrix(0.0)

        # G => diag[-1]   (NxN)
        G = cvxopt.matrix(np.diag(np.ones(numSamples) * -1))

        # h => 0 = [0]    (Nx1)
        h = cvxopt.matrix(np.zeros(numSamples))

        # solves quadratic equation for alpha
        # solution is a dictionary
        solution = cvxopt.solvers.qp(P, q, G, h, A, b)

        # get multipliers, a, from solution
        a = np.ravel(solution['x'])

        # only use lagrange multipliers > 0
        sv = a > 1e-5
        ind = np.arange(len(a))[sv]
        self.a = a[sv]
        self.sv_x = X[sv]   # which rows to keep, ie. which samples
        self.sv_y = y[sv]   # which samples to keep

        # find intercept b from y = sign(w.x - b)
        # not sure where they got this formula from
        self.b = 0
        for n in range(len(self.a)):
            self.b += self.sv_y[n]
            xdotvect = [np.dot(X[ind[n]], X[j]) for j in np.nditer(ind)]
            self.b -= np.sum(self.a * self.sv_y * xdotvect)
        self.b /= len(self.a)

        # find weights w
        self.w = np.zeros(numFeatures)
        for n in range(len(self.a)):
            self.w += self.a[n] * self.sv_y[n] * self.sv_x[n]
    
    
    def predictY(self, X):
        # returns predicted Y for a given X
        # y = sign(w.x + b)
        return np.sign(np.dot(X, self.w) + self.b)