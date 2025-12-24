import numpy as np
import scipy as sp
from scipy.special import *
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'QtAgg', 'GTK3Agg'
import matplotlib.pyplot as plt


# This fcn creates the plots of ce and the eigenvalues
# evidencing the bug in the Mathieu impls.


N = 200                         # Number of sample pts
#v = np.linspace(-180, 180, N)   # Fcn domain -- note this is degrees, not rads.
v = np.linspace(-np.pi, np.pi, N)   # Fcn domain -- this is rads for new impl.


#--------------------------------------------------------
# First make plot of eigenvalues vs. q.  This is an expansion
# of the famous plot DLMF Fig. 28.2.1 to larger m and q values.
# The bug is evidenced by the places where the eigenvalues jump
# from one branch of the plot to another.  This happens when the
# rootfinder jumps to the wrong root.
q = np.linspace(0,50,2500)  # Freq parameter

# Eigenvalues for ce
for m in range(9):
    a = mathieu_a(m,q)
    plt.plot(q, a,'b-')

# Eigenvalues for se
for m in range(1,10):
    b = mathieu_b(m,q)
    plt.plot(q, b,'r--')
plt.title("First Matheiu eigenvalues vs. q")
plt.xlabel("q")
plt.ylabel("Eigenvalues")
    
plt.show()

#--------------------------------------------------------
# Now create a figure and 10 vertically stacked subplots.
# For some values of q we get the wrong Mathieu function.
fig, axs = plt.subplots(10, 1, figsize=(8, 12), sharex=True)

# Array of q (frequency parameter) values.
qs = np.array([25.0, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8, 25.9])
N = len(qs)
m = 6   # Order of Mathieu fcn

for i in range(N):
    q = qs[i]
    ce,_ = mathieu_cem(m,q,v)
    axs[i].plot(v, ce)
    axs[i].grid(True)

axs[-1].set_xlabel('v')
for i in range(N):
    ax = axs[i]
    ax.text(-20,-0.6,'q = %3.1f' % qs[i])        
    
fig.suptitle("Mathieu ce(%d,q,v) for varying q" % m)
plt.show()
