import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import pickle
import logging

# konstanten
a = 4
b = 12
c = 5
E = 21000
N = 5
R = 2.5
r0 = 1.25
order = 1


def save_obj(obj, name):
    with open('obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def radius(x):
    if type(x) == np.ndarray:
        radii = np.zeros(len(x))
        for n, coord in enumerate(x):
            if coord < 0 or coord > a+b+c:
                print("Error, radius(coord) out of bounds!")
            if coord <= a:
                radii[n] = R
            elif coord > a and coord < a+b:
                radii[n] = R + (coord-a)*(r0-R)/b
            else:
                radii[n] = r0
        return radii
    else:
        if x < 0 or x > a+b+c:
            print("Error, radius(x) out of bounds!")
        if x <= a:
            return R
        elif x > a and x < a+b:
            return R + (x-a)*(r0-R)/b
        else:
            return r0


def EA(x):
    return E*np.pi*radius(x)**2


def quad_Gauss(f, a, b, N):
    """
    Approximate the integral of any (callable) function f
    on the interval a,b
    using Gauss Quadrature with N support points.
    Weights/support points are loaded for the reference interval [0,1]
    and must be transformed to the integration interval.
    """
    try:
        weights = load_obj("gauss_lambda")[N]
    except KeyError as e:
        print("Weights for N={} not available".format(N))
        raise e

    try:
        support_points = load_obj("gauss_tau")[N]
    except KeyError as e:
        print("Gausspts for N={} not available".format(N))
        raise e

    logging.debug("quad: Gauss pts: {}".format(support_points))
    # transformation to a-b
    support_points = [a + (b-a)*t for t in support_points]
    function_values = [f(t) for t in support_points]
    logging.debug("quad: Gauss pts transformed: {}".format(support_points))
    logging.debug("quad: f @ Gauss pts: {}".format(function_values))
    logging.debug("quad: weights: {}".format(weights))

    return float((b-a) * np.sum(np.multiply(weights, function_values)))


def phi(x, i, vx):
    """
    Piecewise linear ansatzfunction
    Returns np.ndarray of values phi_i(x)

    x:
    i:  index of ansatzfunktion
    vx: node vector (discretization)
    """

    if not type(x) == np.ndarray:
        x = np.array([x])

    phi = np.zeros(len(x))
    for n, x_ in enumerate(x):
        if i == 0 and x_ == vx[0]:
            phi[n] = 1.0
        elif i == len(vx)-1 and x_ == vx[-1]:
            phi[n] = 1.0
        elif vx[i-1] <= x_ < vx[i]:
            assert i != 0
            phi[n] = (x_ - vx[i-1])/(vx[i] - vx[i-1])
        elif vx[i] <= x_ < vx[i+1]:
            phi[n] = (x_ - vx[i+1])/(vx[i] - vx[i+1])
        else:
            phi[n] = 0
    return phi


def dphi(x, i, vx):
    """
    Derivative of piecewise linear ansatzfunction by x
    Returns np.ndarray of values dphi_i(x)

    x:
    i:  index of ansatzfunktion
    vx: node vector (discretization)
    """

    numpyarray = True
    if not type(x) == np.ndarray:
        x = np.array([x])
        numpyarray = False

    dphi = np.zeros(len(x))
    for n, x_ in enumerate(x):
        if i == 0 and x_ == vx[0]:
            dphi[n] = -1.0
        elif i == len(vx)-1 and x_ == vx[-1]:
            dphi[n] = 1.0
        elif vx[i-1] < x_ <= vx[i]:
            assert i != 0
            dphi[n] = 1/(vx[i] - vx[i-1])
        elif vx[i] < x_ < vx[i+1]:
            dphi[n] = 1/(vx[i] - vx[i+1])
        else:
            dphi[n] = 0

    if not numpyarray:
        return float(dphi)
    else:
        return dphi


def IntegrateFunction(aa, bb, func, order, i, k, vx):
    """
    Integrates func(x)*dphi_i(x)*dphi_k(x) dx on [aa,bb]
    using Gauss Quadrature

    order: specifies number of integration points

    """

    def f(x): return func(x) * dphi(x, i, vx) * dphi(x, k, vx)
    return quad_Gauss(f, aa, bb, N=order)


def NodeStiffness(i, vx, f, order):
    """
    Get contributions to stiffness matrix from node i

    Parameters:
        i: node id
        vx: node vector
        f: function ?
        order: number of integration points to use
    Returns (ii,jj,vv):
            ii: row indices
            jj: column indices
            vv: values
    """
    nvx = len(vx)

    iijjvv = []

    if i == 0:
        a = vx[0]
        b = vx[1]
        v = IntegrateFunction(a, b, f, order, i, 0, vx)
        iijjvv.append((i, 0, v))
        v = IntegrateFunction(a, b, f, order, i, 1, vx)
        iijjvv.append((i, 1, v))
    elif i == nvx-1:
        a = vx[-2]
        b = vx[-1]
        v = IntegrateFunction(a, b, f, order, i, i-1, vx)
        iijjvv.append((i, i-1, v))
        v = IntegrateFunction(a, b, f, order, i, i, vx)
        iijjvv.append((i, i, v))
    else:
        a = vx[i-1]
        b = vx[i]
        v = IntegrateFunction(a, b, f, order, i, i-1, vx)
        iijjvv.append((i, i-1, v))
        v = IntegrateFunction(a, b, f, order, i, i, vx)
        iijjvv.append((i, i, v))

        a = vx[i]
        b = vx[i+1]
        v = IntegrateFunction(a, b, f, order, i, i, vx)
        iijjvv.append((i, i, v))
        v = IntegrateFunction(a, b, f, order, i, i+1, vx)
        iijjvv.append((i, i+1, v))

    for ijv in iijjvv:
        if ijv[2] == 0:
            print(ijv)
            # import ipdb
            # ipdb.set_trace()
    return iijjvv


def BuildStiffnessMatrix(vx, EA, order):
    """
    Gesamtsteifigkeitsmatrix K
    """

    ijv = []
    for node_idx, node_x in enumerate(vx):
        ijv += NodeStiffness(node_idx, vx, EA, order)

    ii, jj, vv = zip(*ijv)

    Kges = scipy.sparse.csr_matrix((vv, (ii, jj)), shape=[len(vx), len(vx)])
    return Kges


def PartitionMatrices(N, idxU):
    """
    Generate paritioning matrics Le, Lf
    N:    total number of nodes
    idxU: nodes indices at which Neumann boundary conditions are defined

    Returns
    Le, Lf matrices so that u = Lf*uf + Le*ue
    """
    Nf = N - len(idxU)
    Ne = len(idxU)

    LfT = scipy.sparse.lil_matrix((Nf, N), dtype='i')
    LeT = scipy.sparse.lil_matrix((Ne, N), dtype='i')

    n_LfT, n_LeT = 0, 0
    for i in range(N):
        if i not in idxU:
            LfT[n_LfT, i] = 1
            n_LfT += 1
        else:
            LeT[n_LeT, i] = 1
            n_LeT += 1

    LfT = LfT.tocsr()
    Lf = LfT.T
    LeT = LeT.tocsr()
    Le = LeT.T

    return Le, Lf


def ReduceSparseMatrix(K, idx):
    """
    Reduces the sparse stiffness matrix K to Kred

    K:   full, sparse stiffness matrix
    idx: node indices at which Neumann boundary conditions are defined

    Returns:
    Kred:  reduced stiffness matrix
    P: Permutation matrix that reduces the node displacement vector accordingly
    """
    N = K.shape[0]
    Le, Lf = PartitionMatrices(N, idx)
    Kred = Lf.T @ K @ Lf

    return Kred


if __name__ == "__main__":
    # plot setup
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)

    # diskretisierung
    n = 0
    vx = np.zeros(N + 2*N + N)  # vektor der knotenkoordinaten
    vx[0: N] = np.arange(0, a, a/N)
    vx[N: N+2*N] = np.arange(a, a+b, b/(2*N))
    vx[N+2*N:] = np.arange(a+b, a+b+c, c/N)
    ax1.scatter(vx, [0]*len(vx), c="orange")

    # darstellung geometrie und knoten
    pltx = np.linspace(0, a+b+c, 100)
    pltgeom = radius(pltx)
    ax1.plot(pltx,  pltgeom, c="blue")
    ax1.plot(pltx, -pltgeom, c="blue")
    # ax1.scatter(vx, [0]*(N+2*N+N), c="orange")

    # gesamtsteifigkeitsmatrix
    Kges = BuildStiffnessMatrix(vx, EA, order)

    # randbedingungen
    # verschiebung
    idxU = [0]
    rbU = [0.0]

    # kraft
    idxF = [len(vx)-1]
    rbF = [1.0]

    # vectoren
    f = scipy.sparse.csc_matrix(
        (rbF, (idxF, [0]*len(idxF))), shape=(len(vx), 1))

    Le, Lf = PartitionMatrices(len(vx), idxU)
    Kred = Lf.T @ Kges @ Lf

    rbU = scipy.sparse.csc_matrix(rbU)
    rhs = Lf.T @ f - Lf.T @ Kges @ Le @ rbU

    # solve
    uf = scipy.sparse.linalg.spsolve(Kred, rhs)
    # weiterverarbeitung als sparse
    uf = scipy.sparse.csc_matrix(uf).T

    # verschiebungen
    u = Le @ rbU + Lf @ uf
    print("Verschiebungen:")
    print(u)

    # kraefte
    fe = Le.T @ Kges @ Le @ rbU + Le.T @ Kges @ Lf @ uf
    f = Le @ fe + Lf @ Lf.T @ f
    print("Kraefte:")
    print(f)

    plt.show()
