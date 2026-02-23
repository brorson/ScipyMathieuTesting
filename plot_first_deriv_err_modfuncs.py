#!/usr/bin/env python
"""
This calls the Mathieu radial fcn impl and gets the fcn value and its
first deriv.  Then it computes the first deriv of the fcn using
a 4th order finite-difference formula.  It then compares the FD result
against the result returned by the fcn impl and makes a heat map of
the error.
"""

import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Convenience fcns.
def modmc1(m, q, x):
    return scipy.special.mathieu_modcem1(m, q, x)
def modmc2(m, q, x):
    return scipy.special.mathieu_modcem2(m, q, x)

def modms1(m, q, x):
    return scipy.special.mathieu_modsem1(m, q, x)
def modms2(m, q, x):
    return scipy.special.mathieu_modsem2(m, q, x)


#==========================================================
# Check first deriv of modmc1
NN = 100
v = np.linspace(1,5,NN)  # Testing domain.
h = 3e-6

tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []
fignum = 1

print("Testing modmc1, h = %e\n" % h)


for q in qs:
    for m in range(MM):  
        #ym3 = modmc1(m, q, v - 3 * h)[0]
        ym2 = modmc1(m, q, v - 2 * h)[0]
        ym1 = modmc1(m, q, v - h)[0]
        y0,y0d  = modmc1(m, q, v)
        yp1 = modmc1(m, q, v + h)[0]
        yp2 = modmc1(m, q, v + 2 * h)[0]
        #yp3 = modmc1(m, q, v + 3 * h)[0]

        # 4th order coeffs for first deriv:
        #  1/12         −2/3    0       2/3     −1/12
        yd = (ym2/12 - 2*ym1/3 + 0*y0 + 2*yp1/3 - yp2/12)/h

        # 6th order coeffs for first deriv:
        #  −1/60        3/20    −3/4    0       3/4     −3/20   1/60
        #yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4 + 0*y0 + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Residual
        r = yd-y0d

        #stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        #l2norm = np.linalg.norm(y0d, ord=2) / np.sqrt(len(y0d))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))

        error.append( e )
        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')
#        print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-15, 5, 5)

plt.figure(fignum)
fignum=fignum+1
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("modmc1 deriv -- rel rms err against 4th order FD")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Check first deriv of modms1
error = []
print("Testing modms1, h = %e\n" % h)

for q in qs:
    for m in range(MM):  
        #ym3 = modms1(m, q, v - 3 * h)[0]
        ym2 = modms1(m, q, v - 2 * h)[0]
        ym1 = modms1(m, q, v - h)[0]
        y0,y0d  = modms1(m, q, v)
        yp1 = modms1(m, q, v + h)[0]
        yp2 = modms1(m, q, v + 2 * h)[0]
        #yp3 = modms1(m, q, v + 3 * h)[0]


        # 4th order coeffs for first deriv:
        #  1/12         −2/3    0       2/3     −1/12
        yd = (ym2/12 - 2*ym1/3 + 0*y0 + 2*yp1/3 - yp2/12)/h

        # 6th order coeffs for first deriv:
        #  −1/60        3/20    −3/4    0       3/4     −3/20   1/60
        #yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4 + 0*y0 + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Residual
        r = yd-y0d

        #stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        #l2norm = np.linalg.norm(y0d, ord=2) / np.sqrt(len(y0d))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))

        error.append( e )

        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')
#        print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)

plt.figure(fignum)
fignum=fignum+1
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("modms1 deriv -- rel rms err against 4th order FD")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Check first deriv of modmc2
v = np.linspace(2,5,NN)  # Testing domain.
error = []
print("Testing modmc2, h = %e\n" % h)

for q in qs:
    for m in range(MM):  
        #ym3 = modmc2(m, q, v - 3 * h)[0]
        ym2 = modmc2(m, q, v - 2 * h)[0]
        ym1 = modmc2(m, q, v - h)[0]
        y0,y0d  = modmc2(m, q, v)
        yp1 = modmc2(m, q, v + h)[0]
        yp2 = modmc2(m, q, v + 2 * h)[0]
        #yp3 = modmc2(m, q, v + 3 * h)[0]

        # 4th order coeffs for first deriv:
        #  1/12         −2/3    0       2/3     −1/12
        yd = (ym2/12 - 2*ym1/3 + 0*y0 + 2*yp1/3 - yp2/12)/h


        # 6th order coeffs for first deriv:
        #  −1/60        3/20    −3/4    0       3/4     −3/20   1/60
        #yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4 + 0*y0 + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Residual
        r = yd-y0d

        #stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        #l2norm = np.linalg.norm(y0d, ord=2) / np.sqrt(len(y0d))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))

        error.append( e )

        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')
#        print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)

plt.figure(fignum)
fignum=fignum+1
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("modmc2 deriv -- rel rms err against 4th order FD")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Check first deriv of modms2
error = []
print("Testing modms2, h = %e\n" % h)

for q in qs:
    for m in range(MM):  
        #ym3 = modms2(m, q, v - 3 * h)[0]
        ym2 = modms2(m, q, v - 2 * h)[0]
        ym1 = modms2(m, q, v - h)[0]
        y0,y0d  = modms2(m, q, v)
        yp1 = modms2(m, q, v + h)[0]
        yp2 = modms2(m, q, v + 2 * h)[0]
        #yp3 = modms2(m, q, v + 3 * h)[0]


        # 4th order coeffs for first deriv:
        #  1/12         −2/3    0       2/3     −1/12
        yd = (ym2/12 - 2*ym1/3 + 0*y0 + 2*yp1/3 - yp2/12)/h

        # 6th order coeffs for first deriv:
        #  −1/60        3/20    −3/4    0       3/4     −3/20   1/60
        #yd = (-ym3/60 + 3*ym2/20 - 3*ym1/4 + 0*y0 + 3*yp1/4 - 3*yp2/20 + yp3/60)/h

        # Residual
        r = yd-y0d

        #stddev = np.std(r)  # Stdev of residual.  Should be zero in ideal case.
        #l2norm = np.linalg.norm(y0d, ord=2) / np.sqrt(len(y0d))
        #e = stddev/l2norm  # Rel err

        # Absolute RMS error
        #e = np.sqrt(np.mean(r**2))
        # Relative RMS error
        e = np.sqrt(np.mean((r/y0)**2))

        error.append( e )

        if not np.isclose( e , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {e}')
#        print(f'm = {m}, q = {q}, error = {e}')

M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)

plt.figure(fignum)
fignum=fignum+1
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis', extend='both')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("modms2 deriv -- rel rms err against 4th order FD")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render




# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
