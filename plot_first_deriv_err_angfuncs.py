#!/usr/bin/env python
"""
This calls the Mathieu angular fcn impl and gets the fcn value and its
first deriv.  Then it computes the first deriv of the fcn using
a 6th order finite-difference formula.  It then compares the FD result
against the result returned by the fcn impl and makes a heat map of
the error.
"""

import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def ce(m, q, x):
    return scipy.special.mathieu_cem(m, q, x)
def se(m, q, x):
    return scipy.special.mathieu_sem(m, q, x)


#==========================================================
# Check first deriv of ce
NN = 100
v = np.linspace(-np.pi,np.pi,NN)
h = 1e-3
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3 = ce(m, q, v - 3 * h)[0]
        ym2 = ce(m, q, v - 2 * h)[0]
        ym1 = ce(m, q, v - h)[0]
        y0,y0d  = ce(m, q, v)
        yp1 = ce(m, q, v + h)[0]
        yp2 = ce(m, q, v + 2 * h)[0]
        yp3 = ce(m, q, v + 3 * h)[0]

 	# 6th order coeffs for first deriv:
        #  −1/60	3/20	−3/4	0	3/4	−3/20	1/60
        yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4 + 0*y0 + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Residual
        r = yd-y0d
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0d, ord=2) / np.sqrt(len(y0d))

        e = stddev/l2norm  # Rel err
        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')
#        print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-20, 10, 5)

plt.figure(1)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Check ce first deriv against finite difference")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Round trip se
#NN = 100
#v = np.linspace(-np.pi,np.pi,NN)
#h = 1e-4
#tol = 1e-3
#MM = 50
#qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3 = se(m, q, v - 3 * h)[0]
        ym2 = se(m, q, v - 2 * h)[0]
        ym1 = se(m, q, v - h)[0]
        y0,y0d  = se(m, q, v)
        yp1 = se(m, q, v + h)[0]
        yp2 = se(m, q, v + 2 * h)[0]
        yp3 = se(m, q, v + 3 * h)[0]

 	# 6th order coeffs:
        #  −1/60	3/20	−3/4	0	3/4	−3/20	1/60
        yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4 + 0*y0 + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Residual
        r = yd-y0d
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0d, ord=2) / np.sqrt(len(y0d))

        e = stddev/l2norm  # Rel err
        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
# levels = np.arange(-20, 20, 5)

plt.figure(2)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Check se first deriv against finite difference")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 

