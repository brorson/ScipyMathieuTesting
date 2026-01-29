#!/usr/bin/env python
"""
This computes the second deriv of the modfcns using a
finite difference formula, and then inserts the second
deriv into the defining ODE and computes the residual
as a test.
"""

import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize



def ce(m, q, x):
    return scipy.special.mathieu_cem(m, q, x)
def se(m, q, x):
    return scipy.special.mathieu_sem(m, q, x)


def mathieu_a(m, q):
    return scipy.special.mathieu_a(m, q)
def mathieu_b(m, q):
    return scipy.special.mathieu_b(m, q)

#==========================================================
# Round trip ce
NN = 100
v = np.linspace(-np.pi,np.pi,NN)
h = 1e-4
h2 = h*h
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3 = ce(m, q, v - 3 * h)[0]
        ym2 = ce(m, q, v - 2 * h)[0]
        ym1 = ce(m, q, v - h)[0]
        y0  = ce(m, q, v)[0] 
        yp1 = ce(m, q, v + h)[0]
        yp2 = ce(m, q, v + 2 * h)[0]
        yp3 = ce(m, q, v + 3 * h)[0]

 	# 6th order coeffs:
        #  1/90 	−3/20 	3/2 	−49/18 	3/2 	−3/20 	1/90
        ydd = (ym3/90 - 3*ym2/20 + 3*ym1/2 - 49*y0/18 + 3*yp1/2 - 3*yp2/20 + yp3/90)/h2

        a = mathieu_a(m,q)

        r = ydd + (a - 2 * q * np.cos(2*v)) * y0
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))

        e = stddev/l2norm  # Rel err
        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-10, 20, 5)

plt.figure(1)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for ce")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Round trip se
NN = 100
v = np.linspace(-np.pi,np.pi,NN)
h = 1e-4
h2 = h*h
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3 = se(m, q, v - 3 * h)[0]
        ym2 = se(m, q, v - 2 * h)[0]
        ym1 = se(m, q, v - h)[0]
        y0  = se(m, q, v)[0] 
        yp1 = se(m, q, v + h)[0]
        yp2 = se(m, q, v + 2 * h)[0]
        yp3 = se(m, q, v + 3 * h)[0]

 	# 6th order coeffs:
        #  1/90 	−3/20 	3/2 	−49/18 	3/2 	−3/20 	1/90
        ydd = (ym3/90 - 3*ym2/20 + 3*ym1/2 - 49*y0/18 + 3*yp1/2 - 3*yp2/20 + yp3/90)/h2

        a = mathieu_b(m,q)

        r = ydd + (a - 2 * q * np.cos(2*v)) * y0
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))

        e = stddev/l2norm  # Rel err
        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-10, 20, 5)

plt.figure(2)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for se")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
