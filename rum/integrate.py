import numpy as np
import scipy.io as sio
import pickle
from numpy.polynomial import polynomial as poly


"""
Functions for numerical integration.

Polynomial notation follows numpy convention:
    Polynomial coefficients in order of increasing degree
    (1, 2, 3) give 1 + 2*x + 3*x**2.

"""


def save_obj(obj, name):
    with open('obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def integrate_polynomial(p, a, b):
    """
    Calculate the exact integral of a polynomial p over bounds (a,b)
    The polynomial is given by a coefficient vector, as in:
    p(x) = 3 + 2x + x³
    coefficients = [3,2,1]
    """
    t_integrals = np.array([(1/(i+1)) * (b**(i+1) - a**(i+1))
                            for i in range(len(p))])
    return sum(p * t_integrals)


def legendre():


if __name__ == "__main__":
    print("Tests")
    print("Test integrate_polynomial...")
    p = np.array([0, 1, 0, 0, -1])
    assert integrate_polynomial(p, 0, 1) == 0.3
