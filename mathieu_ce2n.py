import numpy as np
from scipy.linalg import eigh
from scipy.sparse import diags
from scipy.sparse import diags_array
import matplotlib.pyplot as plt
import scipy

def ce2n(v, q, N, norders):
    for norder in norders:
        j = np.arange(N)
        offsets = [-1,0,1]
        main_diag = np.array((2 * (j)) ** 2)
        main_diag[0] = 0
        off_diag = q * np.ones(N-1)
        off_diag[0] = q
        lower_diag = off_diag.copy()
        lower_diag[0] = 2 * q
        A = diags([lower_diag, main_diag, off_diag], offsets=offsets, shape = (N,N)).toarray()
        A[0,0] = 0
        eig_vals, eig_vecs = eigh(A)
        eig_vecs = eig_vecs / np.linalg.norm(eig_vecs, axis=0)

        coeffs = eig_vecs[:, norder]
        if coeffs[0] < 0:
            coeffs = -coeffs
        C = 2 * coeffs[0] ** 2 + np.sum(coeffs[1:] ** 2)
        coeffs = coeffs / np.sqrt(C)

    return np.array([np.dot(coeffs, np.cos((2 * j) * vk)) for vk in v])

# v here is a list
def se2n(v, q, N, norders): #norders for small n, N is for j 
    for norder in norders:
        j = np.arange(N)
        main_diag = (2 * (j+1)) ** 2
        off_diag = q * np.ones(N-1)
        A = diags([off_diag, main_diag, off_diag], offsets = [-1,0,1], shape = (N, N)).toarray()
        eig_vals, eig_vecs = eigh(A)
        coeffs = eig_vecs[:, norder]
        if coeffs[0] < 0:
            coeffs = -coeffs
        with np.printoptions(precision=13, suppress=True):
            print(eig_vals)
    
    return np.array([np.dot(coeffs, np.sin((2 * j + 2) * vk)) for vk in v])

def sen(v, q, N, norders):
    # check whta is going on
    a=diags_array([q,(2*j+1)**2, q], offsets=[-1,0,1], shape=(N,N)).toarray()
    a = a.copy()
    a[0, 0] = 1 - q
    eval, evec = scipy.linalg.eigh(a)
    evec = evec / np.linalg.norm(evec, axis=0)
    for j in range(N):
        Bo_n = evec[:,j]
        seo=np.zeros_like(v)
        for k in range (n):
            seo += Bo_n[k]*np.sin((2*k+1)*v)
        val = np.argmax(np.abs(seo)>1e-1)
        if seo[val]<0:
            seo=-seo
    # np.set_printoptions(formatter={'float': lambda x: "{0:0.13f}".format(x)})
    # print(eval)
    return seo