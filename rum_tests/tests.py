import unittest
import numpy as np


### WORK IN PROGRESS

from rum.basics import test_spd
from rum.solve import factor_LUP, factor_cholesky, solve_cholesky

class Tests(unittest.TestCase):

    scalar_int = 6
    scalar_float = 2.4
    
    vector = np.array

    A_small = np.array([[0,5,0,4],[-1,-6,0,2],[9,3,-2,1],[-4,7,-1,3]])
    A_spd = np.array([[4.0, 2.0, -1.0, 0.0, 0.0, 0.0],
                    [2.0, 4.0, 1.0, 1.0, 0.0, 1.0],
                    [-1.0, 1.0, 5.0, 3.0, -1.1, 2.0],
                    [0.0, 1.0, 3.0, 4.0, 1.1, 2.5],
                    [0.0, 0.0, -1.1, 1.1, 2.4, 1.0],
                    [0.0, 1.0, 2.0, 2.5, 1.0, 3.0]])

    A_random = np.random.rand(100,100)

    def test_vectornorm(self):
        # Zufallsvektor erstellen
        x = np.random.rand(5)
        p_values = [1,np.inf,2,3.3,4]
        for p in p_values: 
            assert vectornorm(x,p) == np.linalg.norm(x,p)
        
    def test_solve_cholesky(self, ):
        b = np.array([1, -1, 5, 7, 6, -3])
        x = solve_cholesky(A,b)
        assert np.dot(A,x) == b
        
    def test_factor_LUP(self,):
        LU,P = factor_LUP(A)
        
        N = LU.shape[0]

        L = np.identity(N)
        for i in range(0,N):
            for j in range(0,i):
                L[i,j] = LU[i,j]

        U = LU.copy()
        for i in range(0,N):
            for j in range(0,i):
                U[i,j] = 0

        assert np.dot(np.dot(L,U),P) == A

    # direct import
    def test_bookbag(self, ):
        assert Bookbag().mul(1,2,3,4) == 24
