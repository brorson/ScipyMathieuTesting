import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.linalg import eig

def tridiagonal_matrix(n, a, b, c):
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = b[i]
        if i > 0:
            M[i, i-1] = a[i]
        if i < n - 1:
            M[i, i+1] = c[i]
    return M


def make_matrix_e(q, vn):
    N = len(vn)
    h = vn[1]-vn[0]
    c = 2*q*np.cos(2* vn).reshape((len(vn),1))
    hm2 = np.ones((N,1))/(h*h)
    
    
    A = tridiagonal_matrix(N, hm2, -2*hm2-c, hm2)
    A[0,1] = 2*hm2[0]
    A[N-1,N-2] = 2*hm2[0]
    # print(A)
    return A


def mathieu_a(Ne, q, N):
    v = np.linspace(-np.pi, np.pi, N)
    A = make_matrix_e(q, v)
    Draw, k = eig(A)
    Dreal = Draw.real
    idx = np.argsort(Dreal)[::-1]  # Sort in descending order
    
    D = np.zeros(2 * Ne)
    for i in range(Ne):
        # Extract only even order eigenvalues
        D[i] = -Dreal[idx[2*i]]

    return v, D

def mathieu_ce_np(Ne, q, N):
    
    v = np.linspace(-np.pi ,  np.pi, N)
    # A = fourier_diff(N, q)
    A = make_matrix_e(q,v)
   
    Draw, Sraw = np.linalg.eig(A)

    
    Dreal = [val.real for val in Draw]

    Denum = list(enumerate(Dreal))
    Didx = sorted(Denum, key=lambda x: x[1], reverse=True)
    idx = [idx for idx, val in Didx]


    S = Sraw[:, idx]
    ce = S[:, :2*Ne:2]
    
    # Normalisation for the values when v = 0
    zidx = N // 2
    

    for j in range(Ne):
        if ce[zidx, j] < 0:
            ce[:, j] = - ce[:, j]
            

    # Normalisation 28.2.30
    # ce shape : (N, Ne)

    norm_cols = []
    h = v[1] -  v[0]


    # Normalisation for ce_n
    for j in range(ce.shape[1]):
        int_col = 0
        for i in range(1, ce.shape[0]):
            int_f_ij = (0.5) * ((ce[i - 1 , j])**2 + (ce[i, j])**2) * h
            int_col += int_f_ij
        traps_norm = np.sqrt(int_col / np.pi) 
        norm_cols.append(traps_norm)
    
    
    for l in range(ce.shape[1]):
        ce[:, l] = ce[:, l] / norm_cols[l]
    
    return ce 

import numpy as np
import matplotlib.pyplot as plt

def write_mathieu_ce_gvs(q=None):
    """
    This creates a file with golden values in columns used
    to test other impls of the Mathieu ce fcns.
    
    Parameters:
    q : float, optional - if not provided, will prompt for input
    """
    
    if q is None:
        try:
            q = float(input("Enter q value: "))
        except ValueError:
            print("Invalid input for q")
            return
    
    print(f'q = {q:.6f}')
    
    N = 2500  # Number of v values
    # mathieu_ce only operates over the domain [-pi, pi]
    v = np.linspace(-np.pi,  np.pi, N)
    Ne = 35    # Top order of fcn to request.
    
    
    Ss = mathieu_ce_np(Ne, q, N)  # GVs for different orders are
                                # arranged in columns.
   
    
    # Write GVs to a file along with the v value.
    filename = f'/Users/aliamy/Desktop/Mathieu/ScipyMathieuTesting/python_generated_gvs/mathieu_ce_gvs_q{q:.6f}.csv'
    
    with open(filename, 'w') as fh:
        # First write q value to file.
        fh.write(f'{q:14.10f}\n')
        
        # Then write fcn values.
        for i in range(len(v)):
            # Write v value followed by all function values for this v
            line = f'{v[i]}'
            for j in range(Ne):
                line += f', {Ss[i, j].astype(float)}'
            line += '\n'
            fh.write(line)
    
    print(f'Golden values written to {filename}')

# Example usage:
if __name__ == "__main__":
    
    Q =  np.array([-100, -30, -10, -1, -0.1, -0.001, 0.001, 0.01, 0.1, 1, 10, 30, 100])
    # Q =  np.array([ 1])

    for q in Q:
        write_mathieu_ce_gvs(q)

    