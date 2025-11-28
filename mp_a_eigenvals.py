from mpmath import mp, matrix, eig
import matplotlib.pyplot as plt
import numpy as np

mp.dps = 25  # Set precision

"""
This creates the recurrence matrix for the e case and
then finds its eigenvalues and plots them.
"""
 
def tridiagonal_matrix(n, a, b, c):
    M = matrix(n, n)
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
    c = 2*mp.mpf(q)*matrix([mp.cos(mp.mpf(2*x)) for x in vn])
    hm2 = matrix(mp.ones(N,1))/(h*h)
    
    A = tridiagonal_matrix(N, hm2, -2*hm2-c, hm2)
    A[0,1] = 2*hm2[0]
    A[N-1,N-2] = 2*hm2[0]
    return A

def mathieu_a(Ne, q, N):
    v = mp.linspace(-mp.pi, mp.pi, N)
    A = make_matrix_e(q,v)
    Draw, _ = eig(A)
    Dreal = [val.real for val in Draw]
    Denum = list(enumerate(Dreal))
    Didx = sorted(Denum, key=lambda x: x[1], reverse=True)
    idx = [idx for idx, val in Didx]
    D = matrix(2*Ne, 1)
    for i in range(0,Ne):
        # Extract only even order eigenvalues.
        D[i] = -Dreal[idx[2*i]]

    return v, D


#===============================================
if __name__ == "__main__":

    print("Calling mathieu_a ...")
    N = 100     # Length of vn vector
    Nq = 3    # Length of qs vector (number of qs to compute)
    Ne = 6     # Number of eigenvalues to plot
    
    qs = np.linspace(0, 10, Nq)
    D = matrix(Nq,Ne)  
    for i in range(Nq):
        q = qs[i]
        print("q = %f" % q)
        v, AA = mathieu_a(Ne,q,N)
        # Return AA is col vector of eigenvalue orders 0 .. Ne-1
        # D is matrix with eigenvalues for different q arranged
        # in rows
        for j in range(Ne):
            D[i,j] = AA[j,0]
        if (q == mp.mpf(0.0)):
            for j in range(Ne):
                print("a_%d[q=0] = %f" % (j,D[0,j]))


    """
    # Write GVs to a file along with the v value.
    filename = "mathieu_a_gvs.csv"
    fh = open(filename,'w');
    # First write q value to file.
    fh.write('%14.10f\n' % q)
    # Then write fcn values.
    fmt = ['%14.10f, ',repmat('%14.10f, ',[1,Ne-1]),'%f \n'];
    for i=1:length(v):
        fprintf(fh, fmt, v(i), Ss(i,:));
    fh.close()
    """


    # Plot all eigenvalues
    print("Returned.  Now making plot.")
    for i in range(Ne):
        x = [float(qi) for qi in qs]
        y = [float(di) for di in D[:,i]]
        plt.plot(x, y)

    plt.xlabel('q')
    plt.ylabel('eigenvalue a')
    plt.title('Mathieu a eigenvalues')
    plt.show()
