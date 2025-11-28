from mpmath import mp, matrix, eig, sqrt
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser
import os 

# mpmath 
def tridiagonal_matrix(n, a, b, c):
    M = matrix(n, n)
    for i in range(n):
        M[i, i] = b[i]
        if i > 0:
            M[i, i-1] = a[i]
        if i < n - 1:
            M[i, i+1] = c[i]
    return M

def get_middle_elements(c):
    """
    Remove first and last elements from matrix
    Input: (n, 1) matrix
    Output: (n-2, 1) matrix
    
    Equivalent to c[1:-1] in NumPy
    """
    n = c.rows
    if n < 3:
        raise ValueError(f"Matrix must have at least 3 rows, got {n}")
    
    # Create new matrix of size (n-2, 1)
    middle = mp.matrix(n - 2, 1)
    for i in range(n - 2):
        middle[i, 0] = c[i + 1, 0]  # Start from index 1, skip index 0
    
    return middle

# mp math
def make_matrix_o(N, q, vn):

    h = vn[1]-vn[0]
    c = 2*mp.mpf(q)*matrix([mp.cos(mp.mpf(2*x)) for x in vn])
    c = get_middle_elements(c)
    


    hm2 = matrix(mp.ones(N,1))/(h*h)

    
    l = -2*hm2-c


    A = tridiagonal_matrix(N, hm2, l , hm2)


    return A

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

def mpmath_slice_odd_2Ne_explicit(S, Ne):
    """
    Extract odd columns from mpmath matrix
    Equivalent to S[:, 1:2*Ne:2] in NumPy
    Gets columns 1, 3, 5, 7, ... up to 2*Ne-1
    
    Args:
        S: mpmath matrix
        Ne: parameter to determine range (extracts up to column 2*Ne-1)
    
    Returns:
        mpmath matrix with odd columns
    """
    rows = S.rows
    cols_needed = []
    
    # Generate column indices: start=1, stop=2*Ne, step=2
    col = 1  # Start at 1 instead of 0 for odd columns
    while col < 2*Ne and col < S.cols:
        cols_needed.append(col)
        col += 2  # Step by 2 to get odd indices
    
    # Create new matrix
    if not cols_needed:
        return mp.matrix(rows, 0)
    
    co = mp.matrix(rows, len(cols_needed))  # co for "columns odd"
    for i in range(rows):
        for j, col_idx in enumerate(cols_needed):
            co[i, j] = S[i, col_idx]
    
    return co

def mpmath_column_sign_fix(se, zidx, Ne):
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
        diff = se[zidx + 1, j] - se[zidx , j]
        if diff < 0:
            # Flip the sign of entire column j
            for i in range(se.rows):
                se[i, j] = -se[i, j]

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


def add_zero_rows(S):
    """
    Add rows of zeros at top and bottom
    Input: (n, m) matrix
    Output: (n+2, m) matrix
    """
    n, m = S.rows, S.cols
    
    # Create new matrix of size (n+2, m)
    extended = mp.matrix(n + 2, m)
    
    # Copy original matrix to middle rows (rows 1 through n)
    for i in range(n):
        for j in range(m):
            extended[i + 1, j] = S[i, j]
    
    # First and last rows are already zeros by default
    return extended

# mp math


def mathieu_se(Ne, q, N):


    v = mp.linspace(-mp.pi, mp.pi, N)
    A = make_matrix_o(N-2, q,v)

    Draw, Sraw = eig(A) # (N-2, N-2)
   
    Dreal = [val.real for val in Draw]

    Denum = list(enumerate(Dreal))
    Didx = sorted(Denum, key=lambda x: x[1], reverse=True)
    idx = [idx for idx, val in Didx]

    S = extract_columns(Sraw, idx)
    
    se = add_zero_rows(S)

    se = mpmath_slice_odd_2Ne_explicit(se, Ne)
    
    # Getting the correct sign of functions
    zidx = N // 2

    mpmath_column_sign_fix(se, zidx, Ne)
    

    # Normalisation 28.2.30
    # ce shape : (N, Ne)

    h = v[1] -  v[0]

    

    # Normalisation for ce_n
    normalize_mpmath_matrix_explicit(se, h)



    return se # technically some of the normalisation gives me 0.9999999999999999999999999871 see if this is a problem

def write_mathieu_se_gvs(q=None,
                         N=2500,
                         Ne=35,
                         folder='/Users/aliamy/Desktop/Mathieu/ScipyMathieuTesting/python_generated_gvs'
    ):
    """
    This creates a file with golden values in columns used
    to test other impls of the Mathieu se fcns.  
    Parameters:
    q : float, optional - if not provided, will prompt for input
    N : int, Number of v values
    Ne : int, Top order fcn to request
    """
    
    # Read q from user input if not provided
    if q is None:
        try:
            q = float(input("Enter q value: "))
        except ValueError:
            print("Invalid input for q")
            return
    
    print(f'q = {q:.6f}')
    
    # mathieu_ce only operates over the domain [-pi, pi]
    v = mp.linspace(-mp.pi, mp.pi, N)

    # Compute fcn values
    # try:
    Ss = mathieu_se(Ne, q, N)  # GVs for different orders are
                            # arranged in columns.
    
    
    # Write GVs to a file along with the v value.
    filename = f'{folder}/python_mathieu_se_gvs_q{q:.6f}.csv'
    
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
    write_mathieu_se_gvs(q=q, N=N, Ne=Ne, folder=folder)
