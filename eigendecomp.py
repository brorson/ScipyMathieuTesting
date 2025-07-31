from mpmath import mp, matrix, eig

mp.dps = 50

A = matrix([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]])

# Compute eigenvalues and eigenvectors
eigenvals, eigenvecs = eig(A, left=False, right=True)

res = A*eigenvecs[:,0] - eigenvals[0]*eigenvecs[:,0]
print("First eigenvalue residual = \n", res)

res = A*eigenvecs[:,1] - eigenvals[1]*eigenvecs[:,1]
print("Second eigenvalue residual = \n", res)

res = A*eigenvecs[:,2] - eigenvals[2]*eigenvecs[:,2]
print("Third eigenvalue residual = \n", res)


