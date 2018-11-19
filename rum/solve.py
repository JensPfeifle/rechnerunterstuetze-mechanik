from numpy import zeros, sqrt
import numpy as np
import scipy.linalg as la
from .basics import test_spd


def forwardsubs(lower, rhs):
    """
    Solve y in Ly=b with lower triangular matrix L.
    Makes no assumptions about diagonal of L.
    """
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


def solve_thomas(a, b, c, y):
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


def factor_LUP(A):

    eps = 10**-14  # tolerance for zero comparison

    N = A.shape[0]
    LU = A.astype(np.float).copy()
    P = np.identity(N)

    count_pivot = 0

    for j in range(0, N):
        pivot = LU[j, j]
        if abs(pivot) <= 1.0:
            n = np.argmax(abs(LU[j][j:N])) + j  # (shift um i)
            count_pivot += 1
            # swap LU[:, i], LU[:,n]
            LU[:, j], LU[:, n] = LU[:, n], LU[:, j].copy()
            # swap P[i, :], P[m,:]
            P[n, :], P[j, :] = P[j, :], P[n, :].copy()
            pivot = LU[j, j]
        for i in range(j+1, N):
            L = LU[i, j] / pivot
            LU[i] = LU[j]*L - LU[i]
            if abs(LU[i, j]) < eps:
                LU[i, j] = L
            else:
                print("Place for L is not 0! {}".format(LU[i, j]))
                raise ValueError

    return (LU, P)


def solve_LUP(A, b):

    LU, P = factor_LUP(A)
    N = LU.shape[0]

    # LUPx = L(UPx)= Lz
    # solve z
    L = np.identity(N)
    for i in range(0, N):
        for j in range(0, i):
            L[i, j] = LU[i, j]
    z = forwardsubs(L, b)

    # UPx = U(Px) = Uv
    # solve v
    U = LU.copy()
    for i in range(0, N):
        for j in range(0, i):
            U[i, j] = 0
    v = backwardsubs(U, z)

    # v=PX
    x = np.dot(P.T, v)

    return x


def factor_cholesky(A):
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
                sum_ = A[i, j]
                for k in range(j):
                    sum_ = sum_ - L[i, k]*L[j, k]
                if i == j:
                    L[i, j] = sqrt(sum_)
                else:
                    L[i, j] = sum_ / L[j, j]
    return L


def solve_CG(A, y, verbose=True):
    """ Solve Ax=y for a symmetric matric A.

    """

    if not type(A) == np.ndarray:
        print("This function only works with numpy.ndarray types!")
        return None

    EPSILON = 1e-8

    N = A.shape[0]

    x = np.zeros(N)
    p = np.zeros(N)
    d = np.zeros(N)
    r = np.zeros(N)

    KMAX = 1000

    k = 0
    while (k < KMAX):

        if k == 0:  # first iteration
            beta = 0.0
            r = y
        else:
            beta = np.inner(r, r)/rTr_kminus1
        p = r + np.multiply(beta, p)
        d = A.dot(p)
        alpha = np.inner(r, r)/np.inner(d, p)
        x = x + np.multiply(alpha, p)
        rTr_kminus1 = np.inner(r, r)  # altes residium speichern fuer beta
        r = r - np.multiply(alpha, d)

        err = la.norm(r, 2)/la.norm(y, 2)
        k = k + 1

        if (err < EPSILON):
            if verbose:
                print("solve_CG:")
                print("Residual: r = ", la.norm(r, 2))
                print("Took ", k, " iterations.")
            return x

    print("Error. Too many iterations.")
    return None


def solve_cholesky(A, b):
    L = cholesky(A)
    y = forwardsubs(L, b)
    x = backwardsubs(L.T, y)
    return x
