import numpy as np
import scipy as sp
from scipy.special import *
import time

TT = 250


#=======================================================================
# First time computation of eigenvalues
q = np.linspace(0,50,500)  # Freq parameter

start_time = time.time()
for i in range(TT):
    # Eigenvalues for ce
    for m in range(9):
        a = mathieu_a(m,q)
end_time = time.time()
elapsed_time = end_time-start_time
print("TT = %d, Computation of a took %e seconds" % (TT, elapsed_time))

start_time = time.time()
for i in range(TT):
    # Eigenvalues for se
    for m in range(1,10):
        b = mathieu_b(m,q)
end_time = time.time()
elapsed_time = end_time-start_time
print("TT = %d, Computation of b took %e seconds" % (TT, elapsed_time))


#=======================================================================
# Now time computation of function

m =10
q = 25.4
N = 200                         # Number of sample pts
#v = np.linspace(-180, 180, N)   # Fcn domain -- note this is degrees, not rads.
v = np.linspace(-np.pi, np.pi, N)   # Fcn domain -- this is rads for new impl.


start_time = time.time()
for i in range(TT):
    ce,_ = mathieu_cem(m,q,v)
end_time = time.time()
elapsed_time = end_time-start_time
print("TT = %d, Computation of ce took %e seconds" % (TT, elapsed_time))


start_time = time.time()
for i in range(TT):
    se,_ = mathieu_sem(m,q,v)
end_time = time.time()
elapsed_time = end_time-start_time
print("TT = %d, Computation of se took %e seconds" % (TT, elapsed_time))
 
