from numpy import zeros, sqrt
from .basics import test_spd


def thomas(a, b, c, y):
    """
    Solves Ax=y where:

       [[b,c,0,0]
    A=  [a,b,c,0]
        [0,a,b,c]
        [0,0,a,b]]

    a -> lower diagonal elements
    b -> diagonal elements
    c -> upper diagonal elements
    y -> right hand side
    """
    N = len(y)
    u = zeros(N)
    p = zeros(N)

    u[0] = c/b
    p[0] = y[0]/b

    for n in range(1, N):
        u[n] = c/(b-a*u[n-1])
        p[n] = (y[n] - a*p[n-1])/(b-a*u[n-1])

    x = zeros(N)
    x[-1] = p[-1]
    for n in reversed(range(0, N-1)):
        x[n] = p[n] - u[n]*x[n+1]

    return x


def forwardsubs(lower, rhs):
    "Solve y in Ly=b with lower triangular matrix L"
    N = len(rhs)
    y = zeros(N)
    for n in range(N):
        y[n] = (rhs[n] - sum([lower[n, i]*y[i] for i in range(n)]))/lower[n, n]
    return y


def backwardsubs(upper, rhs):
    "Solve x in Rx=y with upper triangular matrix R"
    N = len(rhs)
    x = zeros(N)
    for n in reversed(range(N)):
        x[n] = (rhs[n] - sum([upper[n, i]*x[i]
                              for i in reversed(range(n, N))]))/upper[n, n]
    return x


def cholesky(A, verbose=False):
    """ 
    Cholesky decomposition of A
    Tests if A is spd. 
    Returns lower diagonal matrix L.
    """
    if not test_spd(A):
        return

    L = A.copy()
    N = A.shape[0]
    for i in range(N):
        for j in range(N):
            if j > i:
                L[i, j] = 0
            else:
                if verbose:
                    print("L{}{}=(".format(i, j), end="")
                sum_ = A[i, j]
                if verbose:
                    print("{}".format(A[i, j]), end="")
                for k in range(j):
                    if verbose:
                        print(" - {}*{}".format(L[i, k], L[j, k]),
                              end="")
                    sum_ = sum_ - L[i, k]*L[j, k]
                if i == j:
                    if verbose:
                        print(")->sqrt")
                    L[i, j] = sqrt(sum_)
                else:
                    if verbose:
                        print(")/{}".format(L[j, j]))
                    L[i, j] = sum_ / L[j, j]
    return L


if __name__ == "__main__":
    import numpy as np
    from basics import matprint
    A = np.array([[4.0, 2.0, -1.0, 0.0, 0.0, 0.0],
                  [2.0, 4.0, 1.0, 1.0, 0.0, 1.0],
                  [-1.0, 1.0, 5.0, 3.0, -1.1, 2.0],
                  [0.0, 1.0, 3.0, 4.0, 1.1, 2.5],
                  [0.0, 0.0, -1.1, 1.1, 2.4, 1.0],
                  [0.0, 1.0, 2.0, 2.5, 1.0, 3.0]])

    b = np.array([1, -1, 5, 7, 6, -3])
    L = cholesky(A)
    matprint(L)
    y = forwardsubs(L, b)
    x = backwardsubs(L.T, y)
    print(x)
    print("---")
    print("r={}".format(np.dot(A, x) - b))
