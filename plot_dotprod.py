#!/usr/bin/env python
"""
This computes the dot product of fcns of same type & different order.
"""

import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def ce(m, q, x):
    return scipy.special.mathieu_cem(m, q, x)
def se(m, q, x):
    return scipy.special.mathieu_sem(m, q, x)


# ==========================================================
# Check ce
NN = 1000  # Number of sample pts
#v = np.linspace(-np.pi, np.pi, NN)
v = np.linspace(-180, 180, NN)
h = v[1]-v[0]
# These are failing values of q for ce
#q = 15.55
#q = 17.45
#q = 20.9
q = 25.3
MM = 25 # Highest order to test

error = []
print("Testing ce\n")

for m1 in range(MM):
    y1 = ce(m1, q, v)[0]
    for m2 in range(MM):
        print("m1 = %d, m2 = %d" % (m1, m2))

        y2 = ce(m2, q, v)[0]

        if (m1 == m2):
            e = np.trapz(y1*y2, v)/180   # Should be 1
            #print("e = %e" % e)
        else:
            e = np.trapz(y1*y2, v)     # Should be 0


        error.append(e)

MM1, MM2 = np.meshgrid(range(MM), range(MM))
error_array = np.array(error).reshape(MM1.shape)
levels = np.arange(-15, 5, 5)

eps = 1e-15
log_err = np.log10(np.abs(error_array) + eps)

plt.figure(1)
plt.contourf(MM1, MM2, log_err, levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m1")
plt.ylabel("m2")
plt.title("ce dot product of different orders, q = %f" % q)
plt.draw()  # Draw the plt.figure
plt.pause(0.001)  # Give GUI time to render


# ==========================================================
# Check se
# These are failing values of q for se
q = 31.35
#q = 31.85

error = []
print("Testing se\n")

for m1 in range(MM):
    y1 = se(m1, q, v)[0]
    for m2 in range(MM):
        print("m1 = %d, m2 = %d" % (m1, m2))

        y2 = se(m2, q, v)[0]

        if (m1 == m2):
            e = np.trapz(y1*y2, v)/180   # Should be 1
            #print("e = %e" % e)
        else:
            e = np.trapz(y1*y2, v)     # Should be 0


        error.append(e)

MM1, MM2 = np.meshgrid(range(MM), range(MM))
error_array = np.array(error).reshape(MM1.shape)
levels = np.arange(-15, 5, 5)

log_err = np.log10(np.abs(error_array) + eps)

plt.figure(2)
plt.contourf(MM1, MM2, log_err, levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m1")
plt.ylabel("m2")
plt.title("se dot product of different orders, q = %f" % q)
plt.draw()  # Draw the plt.figure
plt.pause(0.001)  # Give GUI time to render



# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...")
