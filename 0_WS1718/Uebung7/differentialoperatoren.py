
# coding: utf-8

# In[1]:


import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt


# ### Systemmatrix DX zur Ableitung nach x

# In[2]:


def mp_matrix(N):
    h = 1/N
    diag0 = np.repeat([0],N); diag0[0] = -1; diag0[-1] = 1
    diag1 = np.repeat([0.5],N-1); diag1[0] = 1
    diagm1 = np.repeat([-0.5],N-1); diagm1[-1] = -1
    A_mp = sp.diags((diag0, diag1, diagm1),[0,1,-1])
    return A_mp
    return (1/h)*A_mp

def DX(N):
    mp_block = mp_matrix(N)
    DX = sp.block_diag([mp_block.todense()]*N)
    return DX


# ### Systemmatrix DY zur Ableitung nach y

# In[3]:


def DY(N):
    data = np.concatenate([[-1,1],[-0.5,0.5]*(N-2),[-1,1]])
    iidx = np.repeat([N*i for i in range(0,N)],2)
    j_middle = np.concatenate([[N*i,N*(i+2)] for i in range(0,N-2)])
    jidx = np.concatenate([[0,N],j_middle,[N*(N-2),N*(N-1)]])
    
    joined_iidx = np.concatenate([iidx + np.ones(2*N)*n for n in range(N)])
    joined_jidx = np.concatenate([jidx + np.ones(2*N)*n for n in range(N)])
    joined_data = np.tile(data,N)

    DY = sp.coo_matrix((joined_data,(joined_iidx,joined_jidx)))
    return DY


# ### Zweidimensionaler Differenzialoperator DIFF

# In[5]:


def apply_DX(f,N):
    return DX(N).dot(f)

def apply_DY(f,N):
    return DY(N).dot(f)

def DIFF(N):
    return DX(N)+DY(N)

def apply_DIFF(f,N):
    return (DX(N)+DY(N)).dot(f)

def apply_GRAD(f,N):
    return np.array([DX(N).dot(f),DY(N).dot(f)])

def LAPL(N):
    dx = DX(N)
    dy = DY(N)
    return dx.dot(dx) + dy.dot(dy)

def apply_LAPL(f,N):
    dx = DX(N)
    dy = DY(N)
    return (dx.dot(dx) + dy.dot(dy)).dot(f)

def platt(gittermatrix):
    return np.flip(gittermatrix,0).ravel()

def stapel(gittervektor):
    N = int(np.sqrt(gittervektor.size))
    return np.flip(gittervektor.reshape(N,N),0) 
