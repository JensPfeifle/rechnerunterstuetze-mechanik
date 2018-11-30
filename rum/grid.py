import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import copy
from .operators import DX, DY, DIFF, LAPL
from typing import Callable

class UniformGrid():
    
    def __init__(self, xrange, yrange, h):
        self.h = h
        self.N = int(1/h) + 1
        self.xrange = xrange
        self.yrange = yrange
        self.shape = (self.N, self.N)

        self.grid_x = np.linspace(0,xrange,self.N)
        self.grid_y = np.linspace(0,xrange,self.N)
        
    
    def __iter__(self):
        "Unpacks to grid_x, grid_y, h"
        return iter((self.grid_x, self.grid_y, self.h))

    def __getitem__(self,index):
        try:
            i,j = index
        except TypeError as e:
            print("Index must be a tuple of ints, i.e. grid[i,j].")
            raise e
        return (self.grid_x[i], self.grid_y[j])
    
    def plot(self):
        fig, ax = plt.subplots()

        gridpoints = [self[i,j] for i in range(self.N) for j in range(self.N)]
        plot_x = [gp[0] for gp in gridpoints]
        plot_y = [gp[1] for gp in gridpoints]

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        ax.scatter(plot_x, plot_y)
        return plt

import copy



class Field(UniformGrid):

    def __init__(self, xrange, yrange, h,
                 initialvalue=None):
        super().__init__(xrange, yrange, h)

        self.values = np.zeros([self.N, self.N], dtype='d')
        self.boundary = None

        if not type(initialvalue == None):
            if type(initialvalue) == float:
                self.values.fill(float(initialvalue))
            else:
                raise ValueError("initialvalue must be float")

    def asvector(self):
        return self.values.flatten()

    def asarray(self):
        return self.values

    def setboundary(self, indices_i: [], indices_j:[], func: Callable[[float, float],float]):
        """
        Set boundary values of field given i,j indices in the field
        and a function f(x,y) (which could also be constant...)
        """
        assert (len(indices_i) == len(indices_j))
        for i in indices_i:
            for j in indices_j:
                x, y = self.getxy(i, j)
                self.values[i,j] = func(x,y)


    def getxy(self, i, j):
        return (self.grid_x[i], self.grid_y[j])

    def __iter__(self):
        "Unpacks to grid_x, grid_y "
        return iter((self.grid_x, self.grid_y, self.h))

    def __setitem__(self, index, value):
        try:
            i, j = index
        except TypeError as e:
            print("Index must be a tuple of ints, i.e. field[i,j].")
            raise e
        self.values[i, j] = value

    def __getitem__(self, index):
        try:
            i, j = index
        except TypeError as e:
            print("Index must be a tuple of ints, i.e. field[i,j].")
            raise e
        return (self.values[i, j])

    def copy(self):
        return copy.deepcopy(self)

    def plot(self, colormap='cool'):
        fix, ax = plt.subplots()

        gridpoints = [self.getxy(i, j) for i in range(self.N)
                      for j in range(self.N)]
        plot_x = [gp[0] for gp in gridpoints]
        plot_y = [gp[1] for gp in gridpoints]

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        ax.scatter(plot_x, plot_y, c=self.asvector(), cmap=colormap)
        return plt

    def dx(self):
        """
        returns a new field of the same size,
        containing the partial derivative of
        f w.r.t. x
        """
        d_values = DX(self.shape[0]).dot(self.asvector())
        d_field = self.copy()
        d_field.values = d_values.reshape(self.shape)
        return d_field
    
    def dy(self):
        """
        returns a new field of the same size,
        containing the partial derivative of
        f w.r.t. y
        """
        d_values = DY(self.shape[0]).dot(self.asvector())
        d_field = self.copy()
        d_field.values = d_values.reshape(self.shape)
        return d_field
    
    def diff(self):
        """
        returns a new field of the same size,
        containing the derivative of
        f w.r.t. x,y defined as 
        diff(f(x,y) = df(x,y)/dx + df(x,y)/dy
        """
        d_values = DIFF(self.shape[0]).dot(self.asvector())
        d_field = self.copy()
        d_field.values = d_values.reshape(self.shape)
        return d_field
      
    def lapl(self):
        """
        returns a new field of the same size,
        containing the derivative of
        f w.r.t. x,y defined as 
        diff(f(x,y) = df(x,y)/dx + df(x,y)/dy
        """
        d_values = LAPL(self.shape[0]).dot(self.asvector())
        d_field = self.copy()
        d_field.values = d_values.reshape(self.shape)
        return d_field

