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
# Round trip modmc1
NN = 100
u = np.linspace(1,5,NN)
h = 3e-6
h2 = h*h
tol = h # 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        #ym3 = modmc1(m, q, u - 3 * h)[0]
        ym2 = modmc1(m, q, u - 2 * h)[0]
        ym1 = modmc1(m, q, u - h)[0]
        y0  = modmc1(m, q, u)[0] 
        yp1 = modmc1(m, q, u + h)[0]
        yp2 = modmc1(m, q, u + 2 * h)[0]
        #yp3 = modmc1(m, q, u + 3 * h)[0]

        # 4th order coeffs:
        # −1/12 	4/3 	−5/2 	4/3 	−1/12
        ydd = (-ym2/12 + 4*ym1/3 - 5*y0/2 + 4*yp1/3 - yp2/12)/h2

 	# 6th order coeffs:
        #  1/90 	−3/20 	3/2 	−49/18 	3/2 	−3/20 	1/90
        #ydd = (ym3/90 - 3*ym2/20 + 3*ym1/2 - 49*y0/18 + 3*yp1/2 - 3*yp2/20 + yp3/90)/h2

        a = mathieu_a(m,q)

        r = ydd - (a - 2 * q * np.cosh(2*u)) * y0

        #stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        #l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))

        error.append( e )
        
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-15, 5, 5)

plt.figure(1)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modmc1")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render



#==========================================================
# Round trip modmc2
error = []

for q in qs:
    for m in range(MM):    
        #ym3 = modmc2(m, q, u - 3 * h)[0]
        ym2 = modmc2(m, q, u - 2 * h)[0]
        ym1 = modmc2(m, q, u - h)[0]
        y0  = modmc2(m, q, u)[0] 
        yp1 = modmc2(m, q, u + h)[0]
        yp2 = modmc2(m, q, u + 2 * h)[0]
        #yp3 = modmc2(m, q, u + 3 * h)[0]

        # 4th order coeffs:
        # −1/12 	4/3 	−5/2 	4/3 	−1/12
        ydd = (-ym2/12 + 4*ym1/3 - 5*y0/2 + 4*yp1/3 - yp2/12)/h2

 	# 6th order coeffs:
        #  1/90 	−3/20 	3/2 	−49/18 	3/2 	−3/20 	1/90
        #ydd = (ym3/90 - 3*ym2/20 + 3*ym1/2 - 49*y0/18 + 3*yp1/2 - 3*yp2/20 + yp3/90)/h2

        a = mathieu_a(m,q)

        r = ydd - (a - 2 * q * np.cosh(2*u)) * y0

        #stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        #l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))
        
        error.append( e )
        
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)

plt.figure(2)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modmc2")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render

"""
# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 

import sys
sys.exit()
"""

#==========================================================
import math
# Round trip modms1
error = []

for q in qs:
    for m in range(MM):  
        #ym3 = modms1(m, q, u - 3 * h)[0]
        ym2 = modms1(m, q, u - 2 * h)[0]
        ym1 = modms1(m, q, u - h)[0]
        y0  = modms1(m, q, u)[0] 
        yp1 = modms1(m, q, u + h)[0]
        yp2 = modms1(m, q, u + 2 * h)[0]
        #yp3 = modms1(m, q, u + 3 * h)[0]

        # 4th order coeffs:
        # −1/12 	4/3 	−5/2 	4/3 	−1/12
        ydd = (-ym2/12 + 4*ym1/3 - 5*y0/2 + 4*yp1/3 - yp2/12)/h2

 	# 6th order coeffs:
        #  1/90 	−3/20 	3/2 	−49/18 	3/2 	−3/20 	1/90
        #ydd = (ym3/90 - 3*ym2/20 + 3*ym1/2 - 49*y0/18 + 3*yp1/2 - 3*yp2/20 + yp3/90)/h2

        b = mathieu_b(m,q)

        r = ydd - (b - 2 * q * np.cosh(2*u)) * y0

        #stddev = np.std(r)
        #l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))

        error.append( e )

        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)

plt.figure(3)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modms1")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Round trip modms2
error = []

for q in qs:
    for m in range(MM):  
        #ym3 = modms2(m, q, u - 3 * h)[0]
        ym2 = modms2(m, q, u - 2 * h)[0]
        ym1 = modms2(m, q, u - h)[0]
        y0  = modms2(m, q, u)[0] 
        yp1 = modms2(m, q, u + h)[0]
        yp2 = modms2(m, q, u + 2 * h)[0]
        #yp3 = modms2(m, q, u + 3 * h)[0]

        # 4th order coeffs:
        # −1/12 	4/3 	−5/2 	4/3 	−1/12
        ydd = (-ym2/12 + 4*ym1/3 - 5*y0/2 + 4*yp1/3 - yp2/12)/h2

 	# 6th order coeffs:
        #  1/90 	−3/20 	3/2 	−49/18 	3/2 	−3/20 	1/90
        #ydd = (ym3/90 - 3*ym2/20 + 3*ym1/2 - 49*y0/18 + 3*yp1/2 - 3*yp2/20 + yp3/90)/h2

        b = mathieu_b(m,q)

        r = ydd - (b - 2 * q * np.cosh(2*u)) * y0

        #stddev = np.std(r)
        #l2norm = np.linalg.norm(y0, ord=2) / np.sqrt(len(y0))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))
        
        error.append( e )
        
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)

plt.figure(4)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modms2")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
