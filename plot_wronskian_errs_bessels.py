import numpy as np
import scipy as sp
from scipy.special import *

# Use non-GUI backend (no Tk needed)
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================
# This script creates plots of the error of the Bessel
# implementations in SciPy.
# ============================================================


# Parameters to vary
ms = np.arange(1, 51)
MM = len(ms)

# Domain
N = 100
z = np.linspace(0.1, 10, N)


# ============================================================
# Jn & Jn
# ============================================================

print("Computing Wronskian of Jn & Jn")

# Matrices
errs = np.zeros((MM, N))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over m
for i in range(MM):
    m = ms[i]
    print(f"-----------  m = {m}  -----------")

    # Compute Wronskian
    w = jv(m + 1, z) * jv(-m, z) + jv(m, z) * jv(-m - 1, z)
    wtrue = -2 * np.sin(m * np.pi) / (np.pi * z)

    # Error
    err = np.abs(w - wtrue)
    errs[i, :] = np.log10(err + 1.1e-25)

    # Coordinates
    X[i, :] = m
    Y[i, :] = z


# Plot
plt.figure(figsize=(8, 6))

levels = np.arange(-25, 36, 5)

cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)

plt.xlabel("Order m")
plt.ylabel("z")
plt.title("Log10 of Wronskian error -- Jn, Jn")
plt.clim(-20, 10)
plt.colorbar()

# Save figure
plt.savefig("wronskian_Jn_Jn.png", dpi=300, bbox_inches="tight")
plt.close()

print("Saved: wronskian_Jn_Jn.png")


# ============================================================
# Jn & Yn
# ============================================================

print("Computing Wronskian of Jn & Yn")

# Matrices
errs = np.zeros((MM, N))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

# Loop over m
for i in range(MM):
    m = ms[i]
    print(f"-----------  m = {m}  -----------")

    # Compute Wronskian
    w = jv(m + 1, z) * yv(m, z) - jv(m, z) * yv(m + 1, z)
    wtrue = 2 / (np.pi * z)

    # Error
    err = np.abs(w - wtrue)
    errs[i, :] = np.log10(err + 1.1e-25)

    # Coordinates
    X[i, :] = m
    Y[i, :] = z


# Plot
plt.figure(figsize=(8, 6))

levels = np.arange(-25, 36, 5)

cs = plt.contourf(X, Y, errs, levels=levels)
plt.clabel(cs, inline=True)

plt.xlabel("Order m")
plt.ylabel("z")
plt.title("Log10 of Wronskian error -- Jn, Yn")
plt.clim(-20, 10)
plt.colorbar()

# Save figure
plt.savefig("wronskian_Jn_Yn.png", dpi=300, bbox_inches="tight")
plt.close()

print("Saved: wronskian_Jn_Yn.png")


print("\nAll plots generated successfully.")
