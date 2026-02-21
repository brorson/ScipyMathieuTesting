import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_cem, mathieu_a


def fd_deriv(y):
    """
    First derivative using centered finite differences.
    Matches MATLAB-style FD used for verification.
    """
    dydx = np.zeros_like(y)
    dydx[1:-1] = (y[2:] - y[:-2]) / 2.0
    return dydx


def plot_mathieu_ce():
    # ----------------------------------------------------------
    # Plot Mathieu ce functions (Fourier series)
    q = 1.0
    print(f"Plotting Fourier series ce for q = {q}")

    v = np.linspace(0, np.pi / 2, 100)

    # MATLAB default color order
    colors = ['k', 'r', 'g', 'b']  # black, red, green, blue

    # ----------------------------------------------------------
    # Even ce functions: m = 0,2,4,6
    plt.figure(1)
    leg = []

    for i, m in enumerate(range(0, 7, 2)):
        y = mathieu_cem(m, q, v)[0]
        plt.plot(v, y, color=colors[i])
        leg.append(f"m = {m}")

    plt.title("Mathieu ce2m")
    plt.legend(leg, loc="lower left")
    plt.xlabel("v")
    plt.ylabel("ce")

    # ----------------------------------------------------------
    # Even ce function derivs: m = 0,2,4,6
    plt.figure(2)
    leg = []

    for i, m in enumerate(range(0, 7, 2)):
        yd = mathieu_cem(m, q, v)[1]
        plt.plot(v, yd, color=colors[i])
        leg.append(f"m = {m}")

    plt.title("Mathieu ce2m derivative")
    plt.legend(leg, loc="lower left")
    plt.xlabel("v")
    plt.ylabel("ced")
    
    # ----------------------------------------------------------
    # Odd ce functions: m = 1,3,5,7
    plt.figure(3)
    leg = []

    for i, m in enumerate(range(1, 8, 2)):
        y = mathieu_cem(m, q, v)[0]
        plt.plot(v, y, color=colors[i])
        leg.append(f"m = {m}")

    plt.title("Mathieu ce2m+1")
    plt.legend(leg, loc="lower left")
    plt.xlabel("v")
    plt.ylabel("ce")

    # ----------------------------------------------------------
    # Odd ce function derivs: m = 1,3,5,7
    plt.figure(4)
    leg = []

    for i, m in enumerate(range(1, 8, 2)):
        yd = mathieu_cem(m, q, v)[1]
        plt.plot(v, yd, color=colors[i])
        leg.append(f"m = {m}")

    plt.title("Mathieu ce2m+1 derivative")
    plt.legend(leg, loc="lower left")
    plt.xlabel("v")
    plt.ylabel("ced")
    
    # ----------------------------------------------------------
    # Plot ce and derivative together
    N = 250
    u = np.linspace(-np.pi, np.pi, N)
    m = 3
    q = 5.0

    # Compute Mathieu function and derivative
    y, yd = mathieu_cem(m, q, u)

    txt = "m = %d\nq = %5.2f" % (m,q)
    plt.figure(5)
    plt.plot(u, y, label="ce")
    plt.plot(u, yd, label="ced")
    plt.title("ce and ced")
    plt.text(2.5,2.7,txt)
    plt.xlabel("v")
    plt.ylabel("ce, ced")
    plt.legend()

    # ==========================================================
    # Compare against high-order finite-difference approximation
    # DLMF 28.20.11
    m = 7
    N = 1000
    u = np.linspace(0, 2 * np.pi, N)
    h = u[1] - u[0]
    q = 100.0

    # Compute Mathieu function and derivative
    y, yd = mathieu_cem(m, q, u)
    a = mathieu_a(m, q)

    # Second derivative via finite difference
    ydd = fd_deriv(yd) / h

    # Residual of Mathieu equation
    r = ydd + (a - 2 * q * np.cos(2 * u)) * y

    plt_range = slice(5, N - 4)


    # ----------------------------------------------------------
    # Plot deviation from FD approximation
    plt.figure(6)
    plt.plot(u[plt_range], r[plt_range])
    plt.title("Deviation from finite diff approx")
    plt.xlabel("v")
    plt.ylabel("residual")

    # Error diagnostics
    err = np.std(r[plt_range])
    l2norm = np.linalg.norm(r[plt_range])

    print(
        f"Order m = {m}, parameter q = {q}, "
        f"abs err = {err:e}, norm err = {err / l2norm:e}"
    )

    plt.show()


# ----------------------------------------------------------
# Run
if __name__ == "__main__":
    plot_mathieu_ce()
