import numpy as np
import scipy as sp
from scipy.special import *

# Use non-GUI backend (no Tk required)
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================
# This script creates plots of the error of the modified
# Mathieu functions using Wronskian deviation.
# ============================================================


# Parameters
ms = np.arange(1, 51)
qs = np.logspace(-4, 4, 30)

# Domain
N = 100
v = np.linspace(5, 10, N)[:, None]  # Column vector

# True value of Wronskian
wtrue = 2 / np.pi


# ============================================================
# Helper function to plot and save
# ============================================================

def make_plot(fig_num, X, Y, errs, title, filename):

    plt.figure(fig_num, figsize=(8, 6))

    levels = np.arange(-25, 36, 5)

    cs = plt.contourf(X, Y, errs, levels=levels)
    plt.clabel(cs, inline=True)

    plt.xlabel("Order m")
    plt.ylabel("log10(q)")
    plt.title(title)

    plt.clim(-20, 10)
    plt.colorbar()

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {filename}")


# ============================================================
# Mc1 & Mc2
# ============================================================

print("Computing Wronskian of Mc1 & Mc2")

errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

for i, m in enumerate(ms):

    print(f"-----------  m = {m}  -----------")

    for j, q in enumerate(qs):

        y1, y1d = mathieu_modcem1(m, q, v)
        y2, y2d = mathieu_modcem2(m, q, v)

        w = y1 * y2d - y1d * y2

        errs[i, j] = np.log10(np.std((w - wtrue) / wtrue))

        X[i, j] = m
        Y[i, j] = np.log10(q)


make_plot(
    1, X, Y, errs,
    "Log10 of Wronskian rel error -- Mc1 Mc2",
    "wronskian_Mc1_Mc2.png"
)


# ============================================================
# Ms1 & Ms2
# ============================================================

print("Computing Wronskian of Ms1 & Ms2")

errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

for i, m in enumerate(ms):

    print(f"-----------  m = {m}  -----------")

    for j, q in enumerate(qs):

        y1, y1d = mathieu_modsem1(m, q, v)
        y2, y2d = mathieu_modsem2(m, q, v)

        w = y1 * y2d - y1d * y2

        errs[i, j] = np.log10(np.std((w - wtrue) / wtrue))

        X[i, j] = m
        Y[i, j] = np.log10(q)


make_plot(
    2, X, Y, errs,
    "Log10 of Wronskian rel error -- Ms1 Ms2",
    "wronskian_Ms1_Ms2.png"
)


# ============================================================
# Mc1 & Ms2
# ============================================================

print("Computing Wronskian of Mc1 & Ms2")

errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

for i, m in enumerate(ms):

    print(f"-----------  m = {m}  -----------")

    for j, q in enumerate(qs):

        y1, y1d = mathieu_modcem1(m, q, v)
        y2, y2d = mathieu_modsem2(m, q, v)

        w = y1 * y2d - y1d * y2

        errs[i, j] = np.log10(np.std((w - wtrue) / wtrue))

        X[i, j] = m
        Y[i, j] = np.log10(q)


make_plot(
    3, X, Y, errs,
    "Log10 of Wronskian rel error -- Mc1 Ms2",
    "wronskian_Mc1_Ms2.png"
)


# ============================================================
# Ms1 & Mc2
# ============================================================

print("Computing Wronskian of Ms1 & Mc2")

errs = np.zeros((len(ms), len(qs)))
X = np.zeros_like(errs)
Y = np.zeros_like(errs)

for i, m in enumerate(ms):

    print(f"-----------  m = {m}  -----------")

    for j, q in enumerate(qs):

        y1, y1d = mathieu_modsem1(m, q, v)
        y2, y2d = mathieu_modcem2(m, q, v)

        w = y1 * y2d - y1d * y2

        errs[i, j] = np.log10(np.std((w - wtrue) / wtrue))

        X[i, j] = m
        Y[i, j] = np.log10(q)


make_plot(
    4, X, Y, errs,
    "Log10 of Wronskian rel error -- Ms1 Mc2",
    "wronskian_Ms1_Mc2.png"
)


print("\nAll plots generated successfully.")
