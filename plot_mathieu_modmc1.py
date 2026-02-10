import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_modcem1, jv


# ------------------------------------------------------------
# TEMPORARY replacement for mathieu_modmc1
# SciPy does NOT provide mathieu_modmc1
# We use modcem1 so the script runs
# ------------------------------------------------------------
def mathieu_modmc1(m, q, u):
    """
    Placeholder for MATLAB mathieu_modmc1.
    Uses scipy.special.mathieu_modcem1 so code executes.
    Returns only the function values (not derivatives).
    """
    y, _ = mathieu_modcem1(m, q, u)
    return y


def plot_mathieu_modmc1():
    """
    Plots modified Mathieu functions (placeholder implementation).
    """

    # Color order: black, red, green, blue
    colors = ['k', 'r', 'g', 'b']

    # =========================================================
    # Even Mc functions
    q = 1
    v = np.linspace(0, 4, 1000)

    plt.figure(1)
    for i, m in enumerate(range(0, 7, 2)):
        y = mathieu_modmc1(m, q, v)
        plt.plot(v, y, color=colors[i], label=f"m = {m}")

    plt.title("Modified Mathieu of first kind Mc2n (placeholder)")
    plt.xlabel("v")
    plt.ylabel("Function value")
    plt.legend(loc="upper right")
    plt.grid(True)

    # =========================================================
    # Odd Mc functions
    plt.figure(2)
    for i, m in enumerate(range(1, 8, 2)):
        y = mathieu_modmc1(m, q, v)
        plt.plot(v, y, color=colors[i], label=f"m = {m}")

    plt.title("Modified Mathieu of first kind Mc2n+1 (placeholder)")
    plt.xlabel("v")
    plt.ylabel("Function value")
    plt.legend(loc="upper right")
    plt.grid(True)

    # =========================================================
    # Guitarrez paper comparison
    plt.figure(3)
    m = 1
    u = np.linspace(0, 2.5, 100)

    for i, q in enumerate(range(1, 4)):
        # Sign change to match Guitarrez paper
        y = -mathieu_modmc1(m, q, u)
        plt.plot(u, y, color=colors[i], label=f"q = {q}")

    plt.title("modmc1 for varying q (Guitarrez comparison)")
    plt.xlabel("u")
    plt.ylabel("Function value")
    plt.legend()
    plt.grid(True)

    # =========================================================
    # Asymptotic comparison with Bessel function
    # DLMF 28.20.11
    m = 5
    N = 10000
    u = np.linspace(0, 5, N)
    q = 30

    yj = jv(m, 2 * np.sqrt(q) * np.cosh(u))
    ym = mathieu_modmc1(m, q, u)

    plt.figure(4)
    plt.plot(u, yj, label="Bessel J")
    plt.plot(u, ym, label="Mathieu modmc1 (placeholder)")
    plt.title("Asymptotic behavior: Mathieu vs Bessel")
    plt.xlabel("u")
    plt.ylabel("Function value")
    plt.legend()
    plt.grid(True)

    # =========================================================
    # Difference plot
    plt.figure(5)
    diff = ym - yj
    plt.semilogy(u, np.abs(diff))
    plt.title("Difference |Mathieu - Bessel|")
    plt.xlabel("u")
    plt.ylabel("Absolute difference")
    plt.grid(True)

    plt.show()


# =============================================================
# Run the function
# =============================================================
if __name__ == "__main__":
    plot_mathieu_modmc1()