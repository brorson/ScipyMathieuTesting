import numpy as np
import scipy as sp
from scipy.special import *
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'QtAgg', 'GTK3Agg'
import matplotlib.pyplot as plt


# This fcn creates plots of the error of the modified
# Mathieus.  It uses the same computation as used in the
# Matlab impl -- I compute the Wronskian and plot
# its deviation from zero.

# Parameters to vary
ms = np.arange(1, 51)
qs = np.logspace(-4, 4, 30)

# Domain
N = 100
v = np.linspace(5, 10, N)[:, None]  # Column vector

# True value of Wronskian
wtrue = 2/np.pi


#===================================================================
# Mc1 & Mc2
print('Computing Wronskian of Mc1 & Mc2')

# Matrix of error values
errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over q and m
for i in range(len(ms)):
    m = ms[i]
    print(f'-----------  m = {m}  -----------')
    for j in range(len(qs)):
        q = qs[j]
        y1, y1d = mathieu_modcem1(m, q, v)
        y2, y2d = mathieu_modcem2(m, q, v)
        # Compute Wronskian
        w = y1 * y2d - y1d * y2
        
        # Err results
        errs[i, j] = np.log10(np.std(w - wtrue))
        X[i, j] = m
        Y[i, j] = np.log10(q)

plt.figure(1)
levels = np.arange(-25, 36, 5)
cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)
plt.xlabel('Order m')
plt.ylabel('log10(q)')
plt.title('Log10 of Wronskian error -- Mc1 Mc2')
plt.clim(-20, 10)
plt.colorbar()
#plt.show(block=False)
plt.draw()           # Draw the figure
plt.pause(0.001)     # Give GUI time to render

#===================================================================
# Ms1 & Ms2
print('Computing Wronskian of Ms1 & Ms2')
# Matrix of error values
errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over q and m
for i in range(len(ms)):
    m = ms[i]
    print(f'-----------  m = {m}  -----------')
    for j in range(len(qs)):
        q = qs[j]
        y1, y1d = mathieu_modsem1(m, q, v)
        y2, y2d = mathieu_modsem2(m, q, v)
        # Compute Wronskian
        w = y1 * y2d - y1d * y2
        
        # Err results
        errs[i, j] = np.log10(np.std(w - wtrue))
        X[i, j] = m
        Y[i, j] = np.log10(q)

plt.figure(2)
levels = np.arange(-25, 36, 5)
cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)
plt.xlabel('Order m')
plt.ylabel('log10(q)')
plt.title('Log10 of Wronskian error -- ms1 ms2')
plt.clim(-20, 10)
plt.colorbar()
#plt.show(block=False)
plt.draw()           # Draw the figure
plt.pause(0.001)     # Give GUI time to render


#===================================================================
# Mc1 & Ms2
print('Computing Wronskian of Mc1 & Ms2')
# Matrix of error values
errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over q and m
for i in range(len(ms)):
    m = ms[i]
    print(f'-----------  m = {m}  -----------')
    for j in range(len(qs)):
        q = qs[j]
        y1, y1d = mathieu_modcem1(m, q, v)
        y2, y2d = mathieu_modsem2(m, q, v)
        # Compute Wronskian
        w = y1 * y2d - y1d * y2
        
        # Err results
        errs[i, j] = np.log10(np.std(w - wtrue))
        X[i, j] = m
        Y[i, j] = np.log10(q)

plt.figure(3)
levels = np.arange(-25, 36, 5)
cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)
plt.xlabel('Order m')
plt.ylabel('log10(q)')
plt.title('Log10 of Wronskian error -- mc1 ms2')
plt.clim(-20, 10)
plt.colorbar()
#plt.show(block=False)
plt.draw()           # Draw the figure
plt.pause(0.001)     # Give GUI time to render


#===================================================================
# Ms1 & Mc2
print('Computing Wronskian of Ms1 & Mc2')
# Matrix of error values
errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over q and m
for i in range(len(ms)):
    m = ms[i]
    print(f'-----------  m = {m}  -----------')
    for j in range(len(qs)):
        q = qs[j]
        y1, y1d = mathieu_modsem1(m, q, v)
        y2, y2d = mathieu_modcem2(m, q, v)
        # Compute Wronskian
        w = y1 * y2d - y1d * y2
        
        # Err results
        errs[i, j] = np.log10(np.std(w - wtrue))
        X[i, j] = m
        Y[i, j] = np.log10(q)

plt.figure(4)
levels = np.arange(-25, 36, 5)
cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)
plt.xlabel('Order m')
plt.ylabel('log10(q)')
plt.title('Log10 of Wronskian error -- ms1 mc2')
plt.clim(-20, 10)
plt.colorbar()
#plt.show(block=False)
plt.draw()           # Draw the figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
