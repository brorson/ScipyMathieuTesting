#!/usr/bin/env python
"""
This computes the first deriv of the modfcns using a
finite difference formula and compares it to the first
deriv returned from the impl under test.  It then makes
a heat map of the errors.
"""

import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize



def modmc1(m, q, x):
    return scipy.special.mathieu_modcem1(m, q, x)
def modmc2(m, q, x):
    return scipy.special.mathieu_modcem2(m, q, x)

def modms1(m, q, x):
    return scipy.special.mathieu_modsem1(m, q, x)
def modms2(m, q, x):
    return scipy.special.mathieu_modsem2(m, q, x)

def mathieu_a(m, q):
    return scipy.special.mathieu_a(m, q)
def mathieu_b(m, q):
    return scipy.special.mathieu_b(m, q)

#==========================================================
# First deriv test modmc1
NN = 100
v = np.linspace(2,5,NN)
h = 1e-8
#h2 = h*h
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        #ym3        = modmc1(m, q, v - 3 * h)[0]
        ym2        = modmc1(m, q, v - 2 * h)[0]
        ym1        = modmc1(m, q, v - h)[0]
        y0, y0d    = modmc1(m, q, v)
        yp1        = modmc1(m, q, v + h)[0]
        yp2        = modmc1(m, q, v + 2 * h)[0]
        #yp3        = modmc1(m, q, v + 3 * h)[0]

 	# 6th order coeffs for first deriv:
        #  −1/60 	3/20 	−3/4 	0 	3/4 	−3/20 	1/60 	
        #yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4  + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

 	# 4th order coeffs for first deriv:
        #  1/12 	−2/3 	0 	2/3 	−1/12 	
        yd = (ym2/12 - 2*ym1/3  + 2*yp1/3 - yp2/12)/h

        # Compute diff in first derivs
        r = yd -y0d
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))
        e = stddev/l2norm  # Rel err

        #e = np.std(r/yd)

        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')
            """
            if (e > 0.5):
                # Make a plot to see what is wrong.
                plt.figure(99)
                plt.plot(v, y0d, label='y0d')
                plt.plot(v, yd,  label='yd')
                plt.show()
            """
            
M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-10, 20, 5)

plt.figure(1)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("First deriv rel error for modmc1")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 

import sys
sys.exit()

#==========================================================
# First deriv test modms1
NN = 100
v = np.linspace(2,5,NN)
h = 1e-8
#h2 = h*h
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3        = modms1(m, q, v - 3 * h)[0]
        ym2        = modms1(m, q, v - 2 * h)[0]
        ym1        = modms1(m, q, v - h)[0]
        y0, y0d    = modms1(m, q, v)
        yp1        = modms1(m, q, v + h)[0]
        yp2        = modms1(m, q, v + 2 * h)[0]
        yp3        = modms1(m, q, v + 3 * h)[0]

 	# 6th order coeffs for first deriv:
        #  −1/60 	3/20 	−3/4 	0 	3/4 	−3/20 	1/60 	
        yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4  + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Compute diff in first derivs
        r = yd -y0d
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
plt.title("First deriv rel error for modms1")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# First deriv test modmc2
NN = 100
v = np.linspace(2,5,NN)
h = 1e-8
#h2 = h*h
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3        = modmc2(m, q, v - 3 * h)[0]
        ym2        = modmc2(m, q, v - 2 * h)[0]
        ym1        = modmc2(m, q, v - h)[0]
        y0, y0d    = modmc2(m, q, v)
        yp1        = modmc2(m, q, v + h)[0]
        yp2        = modmc2(m, q, v + 2 * h)[0]
        yp3        = modmc2(m, q, v + 3 * h)[0]

 	# 6th order coeffs for first deriv:
        #  −1/60 	3/20 	−3/4 	0 	3/4 	−3/20 	1/60 	
        yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4  + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Compute diff in first derivs
        r = yd -y0d
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))

        e = stddev/l2norm  # Rel err
        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-10, 20, 5)

plt.figure(3)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("First deriv rel error for modmc2")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# First deriv test modms2
NN = 100
v = np.linspace(2,5,NN)
h = 1e-8
#h2 = h*h
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        ym3        = modms2(m, q, v - 3 * h)[0]
        ym2        = modms2(m, q, v - 2 * h)[0]
        ym1        = modms2(m, q, v - h)[0]
        y0, y0d    = modms2(m, q, v)
        yp1        = modms2(m, q, v + h)[0]
        yp2        = modms2(m, q, v + 2 * h)[0]
        yp3        = modms2(m, q, v + 3 * h)[0]

 	# 6th order coeffs for first deriv:
        #  −1/60 	3/20 	−3/4 	0 	3/4 	−3/20 	1/60 	
        yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4  + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Compute diff in first derivs
        r = yd -y0d
        stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))

        e = stddev/l2norm  # Rel err
        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-10, 20, 5)

plt.figure(4)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("First deriv rel error for modms2")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
