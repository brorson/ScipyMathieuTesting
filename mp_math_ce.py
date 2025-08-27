from mpmath import mp, matrix, eig, sqrt,mpf
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser
import os 

def tridiagonal_matrix(n, a, b, c):
    '''
    a: upper major diagonal
    b: major diagonal
    c: lower major diagonal
    '''

    M = matrix(n, n)

    # print(f'shape a = {a.shape}, shape b = {b.shape}, c shape = {c.shape}')
    for i in range(n):
        M[i, i] = b[i]
        if i > 0:
            M[i, i-1] = a[i]
        if i < n - 1:
            M[i, i+1] = c[i]

    return M


def make_matrix_e(q, vn):
    '''
    q: q parameter for Mathieu Function
    vn: input values
    '''
    N = len(vn)
    h = vn[1]-vn[0]
    c = 2*mp.mpf(q)*matrix([mp.cos(mp.mpf(2*x)) for x in vn])
    hm2 = matrix(mp.ones(N, 1))/(h*h)
    A = tridiagonal_matrix(N, hm2, -2*hm2-c, hm2)
    A[0, 1] = 2*hm2[0]
    A[N-1, N-2] = 2*hm2[0]
    return A


def mathieu_a(Ne, q, N):
    v = mp.linspace(-np.pi, np.pi, N)
    A = make_matrix_e(q, v)
    Draw, _ = eig(A)
    # print(len(Draw))
    Dreal = [val.real for val in Draw]
    Denum = list(enumerate(Dreal))
    Didx = sorted(Denum, key=lambda x: x[1], reverse=True)
    idx = [idx for idx, val in Didx]
    D = matrix(2*Ne, 1)
    for i in range(0, Ne):
        # Extract only even order eigenvalues.
        D[i] = -Dreal[idx[2*i]]

    return v, D

# -----------------------
#   MP Math utilities
# -----------------------


def extract_single_column(Sraw, col_idx):
    """Extract a single column from mpmath matrix"""
    rows = Sraw.rows
    S = mp.matrix(rows, 1)
    for i in range(rows):
        S[i, 0] = Sraw[i, col_idx]
    return S


def extract_columns(Sraw, idx):
    """Extract multiple columns from mpmath matrix using list of indices"""
    rows = Sraw.rows
    if isinstance(idx, int):
        # Single column case
        return extract_single_column(Sraw, idx)    
    # Multiple columns case
    cols = len(idx)
    S = mp.matrix(rows, cols)
    for i in range(rows):
        for j, col_idx in enumerate(idx):
            S[i, j] = Sraw[i, col_idx]
    return S


def mpmath_slice_even_2Ne_explicit(S, Ne):
    """Extract even columns until index 2Ne"""
    rows = S.rows
    cols_needed = []
    
    # Generate column indices: start=0, stop=2*Ne, step=2
    col = 0
    while col < 2*Ne and col < S.cols:
        cols_needed.append(col)
        col += 2
    # Create new matrix
    if not cols_needed:
        return mp.matrix(rows, 0)  
    ce = mp.matrix(rows, len(cols_needed))
    for i in range(rows):
        for j, col_idx in enumerate(cols_needed):
            ce[i, j] = S[i, col_idx]
    return ce


def mpmath_column_sign_fix(ce, zidx, Ne):
    """
    mpmath equivalent of:
    for j in range(Ne):
        if ce[zidx, j] < 0:
            ce[:, j] = - ce[:, j]
    Args:
        ce: mpmath matrix
        zidx: row index to check for sign
        Ne: number of columns to process
    """
    for j in range(Ne):
        if ce[zidx, j] < 0:
            # Flip the sign of entire column j
            for i in range(ce.rows):
                ce[i, j] = -ce[i, j]


def normalize_mpmath_matrix_explicit(ce, h):
    """
    More explicit version ensuring all operations use mpmath
    """
    norm_cols = []  
    # Convert h to mpmath if it isn't already
    h_mp = mp.mpf(h)  
    # First loop: calculate normalization factors
    for j in range(ce.cols):
        int_col = mp.mpf(0)
        for i in range(1, ce.rows):
            # Ensure all operations are in mpmath
            term1 = ce[i - 1, j] ** 2
            term2 = ce[i, j] ** 2
            int_f_ij = mp.mpf(0.5) * (term1 + term2) * h_mp
            int_col += int_f_ij
        traps_norm = mp.sqrt(int_col / mp.pi) 
        norm_cols.append(traps_norm) 
    # Second loop: normalize each column
    for l in range(ce.cols):
        norm_factor = norm_cols[l]
        for i in range(ce.rows):
            ce[i, l] = ce[i, l] / norm_factor


def mathieu_ce(Ne, q, N):
    v = mp.linspace(-mp.pi, mp.pi, N)
    A = make_matrix_e(q, v)

    Draw, Sraw = eig(A)
    Dreal = [val.real for val in Draw]

    Denum = list(enumerate(Dreal))
    Didx = sorted(Denum, key=lambda x: x[1], reverse=True)
    idx = [idx for idx, val in Didx]

    S = extract_columns(Sraw, idx)
    ce = mpmath_slice_even_2Ne_explicit(S, Ne)

    # Normalisation for the values when v = 0
    zidx = N // 2
    mpmath_column_sign_fix(ce, zidx, Ne)

    # Normalisation 28.2.30
    # ce shape : (N, Ne)

    h = v[1] - v[0]

    normalize_mpmath_matrix_explicit(ce, h)
    print(f'ce mp = {ce}')
    return ce


#Final function to write the ce golden values
def write_mathieu_ce_gvs(q=None,
                         N=2500,
                         Ne=35,
                         folder='/Users/aliamy/Desktop/Mathieu/ScipyMathieuTesting/python_generated_gvs'
    ):
    """
    This creates a file with golden values in columns used
    to test other impls of the Mathieu ce fcns.  
    Parameters:
    q : float, optional - if not provided, will prompt for input
    N : int, Number of v values
    Ne : int, Top order fcn to request
    """
    if q is None:
        try:
            q = float(input("Enter q value: "))
        except ValueError:
            print("Invalid input for q")
            return
    print(f'q = {q:.6f}')
    # mathieu_ce only operates over the domain [-pi, pi]
    v = mp.linspace(-mp.pi, mp.pi, N)
    Ss = mathieu_ce(Ne, q, N)  # GVs for different orders are arranged in columns.
    # Write GVs to a file along with the v value.
    filename = f'{folder}/python_mathieu_ce_gvs_q{q:.6f}.csv'   
    with open(filename, 'w') as fh:
        # First write q value to file.
        fh.write(f'{q:14.10f}\n')
        # Then write fcn values.
        for i in range(len(v)):
            # Write v value followed by all function values for this v
            line = f'{v[i]}'
            for j in range(Ne):
                line += f', {Ss[i, j]}'
            line += '\n'
            fh.write(line)
    print(f'Golden values written to {filename}')


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-q', default=0.0, type=float)
    parser.add_argument('-N', default=2500, type=int)
    parser.add_argument('-Ne', default=35, type=int)
    parser.add_argument('--folder', type=os.path.abspath)
    parser.add_argument('--dps', type=int, default=25)
    args = parser.parse_args()

    mp.dps = 25  # Set precision
    
    q, N, Ne = args.q, args.N, args.Ne
    folder = args.folder
    write_mathieu_ce_gvs(q=q, N=N, Ne=Ne, folder=folder)
