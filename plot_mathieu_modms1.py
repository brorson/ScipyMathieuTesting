import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_modsem1, jv


# ------------------------------------------------------------
# TEMPORARY replacement for mathieu_modms1
# SciPy does NOT provide mathieu_modms1
# Uses modsem1 as a placeholder so code runs
# ------------------------------------------------------------
def mathieu_modms1(m, q, u):
    """
    Placeholder for MATLAB mathieu_modms1.
    Uses scipy.special.mathieu_modsem1 so the script executes.
    Returns only the function values (not derivatives).
    """
    y, _ = mathieu_modsem1(m, q, u)
    return y


def plot_mathieu_modms1():
    """
    Plots modified Mathieu Ms functions of the first kind
    (placeholder implementation).
    """
    # Color order: black, red, green, blue
    colors = ['k', 'r', 'g', 'b']

    # =========================================================
    # Even Ms functions (m = 2,4,6,8)
    q = 1
    v = np.linspace(0, 4, 1000)

    plt.figure(1)
    for i, m in enumerate(range(2, 9, 2)):
        y = mathieu_modms1(m, q, v)
        plt.plot(v, y, color=colors[i], label=f"m = {m}")

    plt.title("Modified Mathieu of first kind Ms2n (placeholder)")
    plt.xlabel("v")
    plt.ylabel("Function value")
    plt.legend()
    plt.grid(True)

    # =========================================================
    # Odd Ms functions (m = 1,3,5,7)
    plt.figure(2)
    for i, m in enumerate(range(1, 8, 2)):
        y = mathieu_modms1(m, q, v)
        plt.plot(v, y, color=colors[i], label=f"m = {m}")

    plt.title("Modified Mathieu of first kind Ms2n+1 (placeholder)")
    plt.xlabel("v")
    plt.ylabel("Function value")
    plt.legend()
    plt.grid(True)

    # =========================================================
    # Guitarrez paper comparison
    plt.figure(3)
    m = 1
    u = np.linspace(0, 2.5, 100)

    for i, q in enumerate(range(1, 4)):
        y = mathieu_modms1(m, q, u)
        plt.plot(u, y, color=colors[i], label=f"q = {q}")

    plt.title("modms1 for varying q (Guitarrez comparison)")
    plt.xlabel("u")
    plt.ylabel("Function value")
    plt.legend()
    plt.grid(True)

    # =========================================================
    # Asymptotic comparison with Bessel J
    # DLMF 28.20.11
    plt.figure(4)
    m = 5
    N = 10000
    u = np.linspace(0, 5, N)
    q = 1

    yj = jv(m, 2 * np.sqrt(q) * np.cosh(u))
    ym = mathieu_modms1(m, q, u)

    plt.plot(u, yj, label="Bessel J")
    plt.plot(u, ym, label="Mathieu modms1 (placeholder)")
    plt.title("Asymptotic behavior: modms1 vs Bessel J")
    plt.xlabel("u")
    plt.ylabel("Function value")
    plt.legend()
    plt.grid(True)

    # =========================================================
    # Difference plot
    plt.figure(5)
    diff = ym - yj
    plt.semilogy(u, np.abs(diff))
    plt.title("Difference |Mathieu - Bessel J|")
    plt.xlabel("u")
    plt.ylabel("Absolute difference")
    plt.grid(True)

    plt.show()


# =============================================================
# Run the function
# =============================================================
if __name__ == "__main__":
    plot_mathieu_modms1()