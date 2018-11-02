from numpy import trace
from numpy.linalg import det, matrix_rank, eig, eigh


def matprops(A):
    print("detA = " + str(det(A)))
    print("rgA = " + str(matrix_rank(A)))
    print("symA = \n" + str(1/2*(A + A.T)))
    print("skwA = \n" + str(1/2*(A - A.T)))
    print("trA = " + str(trace(A)))
    print("lambda_i:\n{} \n u_i: \n{}".format(eig(A)[0], eig(A)[1]))


def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x))
                      for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            fmtstr = "{:"+str(col_maxes[i])+fmt+"}"
            print(fmtstr.format(y), end=" ")
        print("")


def vectornorm(v, p):
    d = np.size(v)

    if p == 1:
        return sum(v)

    elif p == np.inf:
        return max(v)

    elif p == 2:
        return (sum(x**2))**(1/2)

    else:
        n = 0
        for i in range(0, d):
            n = n + abs(v[i])**p
        return n**(1/p)


def test_spd(A):
    """ 
    Returns True if the A is symmetric positive-definite (spd)
    else returns False
    """
    # symmetry
    if not (A.transpose() == A).all():
        return False
    # positive eigenvalues
    lamdas, _ = eigh(A)
    for l in lamdas:
        if l <= 0:
            return False
    return True


if __name__ == '__main__':
    import numpy as np
    from numpy.linalg import norm
    print("Test vectornorm:")
    # Zufallsvektor erstellen
    x = np.random.rand(5)
    p_values = [1,np.inf,2,3.3,4]
    print("Zufallsvektor: {}".format(x))
    for p in p_values:
        print("p={}".format(p))      
        print("my norm: {}".format(vectornorm(x,p)))
        print("la.norm: {}".format(norm(x,p)))
