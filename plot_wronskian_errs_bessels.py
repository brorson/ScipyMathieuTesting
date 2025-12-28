import numpy as np
import scipy as sp
from scipy.special import *
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'QtAgg', 'GTK3Agg'
import matplotlib.pyplot as plt

# This fcn creates plots of the error of the Bessel impls in
# scipy.  Scipy seems to use the Cepehes library for the Bessels.

# Parameters to vary
ms = np.arange(1, 51)
MM = len(ms)

# Domain
N = 100
z = np.linspace(.1, 10, N)  # Column vector


#===================================================================
# Jn & Jn
print('Computing Wronskian of Jn & Jn')

# Matrix of error values
errs = np.zeros((MM, N))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over m
for i in range(len(ms)):
    m = ms[i]
    print(f'-----------  m = {m}  -----------')

    # Compute Wronskian -- this is mth column
    w = jv(m+1,z)*jv(-m,z) + jv(m,z)*jv(-m-1,z)
    wtrue = -2*np.sin(m*np.pi)/(np.pi*z)

    # Compute pointwise error at each z
    err = np.abs(w - wtrue)
    errs[i, :] = np.log10(err + 1.1e-25)  # Add small value to avoid log(0)
    
    # Fill coordinate arrays
    X[i, :] = m
    Y[i, :] = z


plt.figure(1)
levels = np.arange(-25, 36, 5)
cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)
plt.xlabel('Order m')
plt.ylabel('z')
plt.title('Log10 of Wronskian error -- Jn, Jn')
plt.clim(-20, 10)
plt.colorbar()
#plt.show(block=False)
plt.draw()           # Draw the figure
plt.pause(0.001)     # Give GUI time to render


#===================================================================
# Jn & Yn
print('Computing Wronskian of Jn & Yn')

# Matrix of error values
errs = np.zeros((MM, N))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over m
for i in range(len(ms)):
    m = ms[i]
    print(f'-----------  m = {m}  -----------')

    # Compute Wronskian -- this is mth column
    w = jv(m+1,z)*yv(m,z) - jv(m,z)*yv(m+1,z)
    wtrue = 2/(np.pi*z)

    # Compute pointwise error at each z
    err = np.abs(w - wtrue)
    errs[i, :] = np.log10(err + 1.1e-25)  # Add small value to avoid log(0)
    
    # Fill coordinate arrays
    X[i, :] = m
    Y[i, :] = z


plt.figure(2)
levels = np.arange(-25, 36, 5)
cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)
plt.xlabel('Order m')
plt.ylabel('z')
plt.title('Log10 of Wronskian error -- Jn, Yn')
plt.clim(-20, 10)
plt.colorbar()
#plt.show(block=False)
plt.draw()           # Draw the figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
