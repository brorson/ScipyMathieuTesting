import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_modsem2, yv


# ------------------------------------------------------------
# TEMPORARY replacement for mathieu_modms2
# SciPy does NOT provide mathieu_modms2
# Uses modsem2 as a placeholder so code runs
# ------------------------------------------------------------
def mathieu_modms2(m, q, u):
    """
    Placeholder for MATLAB mathieu_modms2.
    Uses scipy.special.mathieu_modsem2 so the script executes.
    Returns only the function values (not derivatives).
    """
    y, _ = mathieu_modsem2(m, q, u)
    return y


def plot_mathieu_modms2():
    """
    Plots modified Mathieu Ms functions of the second kind
    (placeholder implementation).
    """
    # Color order: black, red, green, blue
    colors = ['k', 'r', 'g', 'b']

    q = 1
    vmax = 4
    v = np.linspace(0, vmax, 1000)

    # =========================================================
    # Even Ms functions (m = 2,4,6,8)
    plt.figure(1)
    for i, m in enumerate(range(2, 9, 2)):
        y = mathieu_modms2(m, q, v)
        plt.plot(v, y, color=colors[i], label=f"m = {m}")

    plt.title("Modified Mathieu of second kind Ms2n")
    plt.legend()
    plt.ylim([-2, 2])

    # =========================================================
    # Odd Ms functions (m = 1,3,5,7)
    plt.figure(2)
    for i, m in enumerate(range(1, 8, 2)):
        y = mathieu_modms2(m, q, v)
        plt.plot(v, y, color=colors[i], label=f"m = {m}")

    plt.title("Modified Mathieu of second kind Ms2n+1")
    plt.legend()
    plt.ylim([-2, 2])

    # =========================================================
    # Asymptotic comparison with Bessel Y
    # DLMF 28.20.11
    plt.figure(4)
    m = 7
    N = 1000
    u = np.linspace(0, vmax, N)
    q = 1

    yy = yv(m, 2 * np.sqrt(q) * np.cosh(u))
    ym = mathieu_modms2(m, q, u)

    # Sign change to match Bessel function (commented out in MATLAB)
    # if np.sign(ym[-1]) != np.sign(yy[-1]):
    #     ym = -ym

    plt.plot(u, yy, label="Bessel Y")
    plt.plot(u, ym, label="Mathieu modms2")
    plt.title("Asymptotic behavior: modms2 compared to Bessel Y")
    plt.legend()
    plt.ylim([-6, 6])

    # =========================================================
    # Difference plot
    plt.figure(5)
    diff = ym - yy
    plt.semilogy(u, np.abs(diff))
    plt.title("Difference")

    plt.show()


# =============================================================
# Run the function
# =============================================================
if __name__ == "__main__":
    plot_mathieu_modms2()