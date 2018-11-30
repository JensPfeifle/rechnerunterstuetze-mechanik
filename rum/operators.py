import numpy as np
from scipy import sparse


def DX(N):
    diag0 = np.repeat([0], N)
    diag0[0] = -1
    diag0[-1] = 1
    diag1 = np.repeat([0.5], N-1)
    diag1[0] = 1
    diagm1 = np.repeat([-0.5], N-1)
    diagm1[-1] = -1
    mp_block = sparse.diags((diag0, diag1, diagm1), [0, 1, -1])
    DX = sparse.block_diag([mp_block.todense()]*N)
    return DX


def DY(N):
    data = np.concatenate([[-1, 1], [-0.5, 0.5]*(N-2), [-1, 1]])
    iidx = np.repeat([N*i for i in range(0, N)], 2)
    j_middle = np.concatenate([[N*i, N*(i+2)] for i in range(0, N-2)])
    jidx = np.concatenate([[0, N], j_middle, [N*(N-2), N*(N-1)]])

    joined_iidx = np.concatenate(
        [iidx + np.ones(2*N, dtype=np.int_)*n for n in range(N)])
    joined_jidx = np.concatenate(
        [jidx + np.ones(2*N, dtype=np.int_)*n for n in range(N)])
    joined_data = np.tile(data, N)
    DY = sparse.coo_matrix((joined_data, (joined_iidx, joined_jidx)))
    return -1*DY


def DIFF(N):
    return DX(N)+DY(N)

def LAPL(N):
    dx = DX(N)
    dy = DY(N)
    return dx.dot(dx) + dy.dot(dy)