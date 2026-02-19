import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


def tridiagonal_matrix_np(n, a, b, c):
    M = np.zeros((n, n))
    # print(f'shape a = {a.shape}, shape b = {b.shape}, c shape = {c.shape}')
    for i in range(n):
        M[i][i] = b[i]
        if i > 0:
            M[i][i-1] = a[i]
        if i < n - 1:
            M[i][i+1] = c[i]

    return M

def make_matrix_o_np(N, q, vn):
    
    h = vn[1]-vn[0]
    c = 2*q*np.cos(2* vn[1:-1]).reshape(len(vn)- 2,1)
    
    
    hm2 = np.ones((N,1))/(h*h)
    l = -2*hm2-c
   
    
    A = tridiagonal_matrix_np(N, hm2, l , hm2)

    # No boundary condition
   
    # print(A)
    return A


def mathieu_se_np(Ne, q, N):


    v = np.linspace(0, 2 * np.pi, N)
    A = make_matrix_o_np(N-2, q,v)

    # A = fourier_diff(N, q)


    # print(f'A= {A}')
    Draw, Sraw = np.linalg.eig(A) # (N-2, N-2)
    # print(f'Sraw = {Sraw}')
    # print(f'D = {Draw}')
    Dreal = [val.real for val in Draw]

    Denum = list(enumerate(Dreal))
    Didx = sorted(Denum, key=lambda x: x[1], reverse=True)
    idx = [idx for idx, val in Didx]


    S = Sraw[:, idx]
    
    

    se = np.vstack([np.zeros((1, S.shape[1])), S, np.zeros((1, S.shape[1]))])

    se = se[:, 1:2*Ne:2]
    
    # Getting the correct sign of functions
    zidx = N // 2
    # h = v[1] - v[0]


    for j in range(Ne):
        diff = se[zidx + 1, j] - se[zidx , j]
        if diff < 0:
            se[:, j] = - se[:, j]
            

    # Normalisation 28.2.30
    # ce shape : (N, Ne)

    norm_cols = []
    h = v[1] -  v[0]

   
    # Normalisation for ce_n
    for j in range(se.shape[1]):
        int_col = 0
        for i in range(1, se.shape[0]):
            int_f_ij = (0.5) * ((se[i - 1 , j])**2 + (se[i, j])**2) * h
            int_col += int_f_ij
        # int_col += (0.5) * ((ce[i - 1 , j])**2 + (ce[i, j])**2) * h
        traps_norm = np.sqrt(int_col / np.pi) 
        norm_cols.append(traps_norm)
    
    
    for l in range(se.shape[1]):
        # print(norm_cols[l])
        se[:, l] = se[:, l] / norm_cols[l]
    
    print(f'se after norm = {se}')



    return se 




def write_mathieu_se_gvs(q=None):
    """
    This creates a file with golden values in columns used
    to test other impls of the Mathieu ce fcns.
    
    Parameters:
    q : float, optional - if not provided, will prompt for input
    """
    
    # Read q from user input if not provided
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
    
    # Compute fcn values
    # try:
    Ss = mathieu_se_np(Ne, q, N)  # GVs for different orders are
                                # arranged in columns.
    
    
    # Write GVs to a file along with the v value.
    filename = f'/Users/aliamy/Desktop/Mathieu/ScipyMathieuTesting/python_generated_gvs/mathieu_se_gvs_q{q:.6f}.csv'
    
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

# Example usage:
if __name__ == "__main__":
    
    
    # Or without arguments (will prompt for input):
    Q =  np.array([-100, -30, -10, -1, -0.1, -0.001, 0.001, 0.01, 0.1, 1, 10, 30, 100])
    # Q =  np.array([0.1])

    for q in Q:
        write_mathieu_se_gvs(q)

    # incorrect value written out
    