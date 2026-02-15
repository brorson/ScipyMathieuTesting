import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_sem


def plot_mathieu_se():
    # ----------------------------------------------------------
    # This plots the Mathieu se functions

    q = 1.0
    v = np.linspace(0, np.pi / 2, 100)

    # Color order: black, red, green, blue
    colors = ['k', 'r', 'g', 'b']

    # ----------------------------------------------------------
    # First: even se functions (m = 2,4,6,8)
    plt.figure(1)
    leg = []

    for i, m in enumerate(range(2, 9, 2)):
        # MATLAB: y = mathieu_se(m, q, v)
        # SciPy:  mathieu_sem(m, q, v)[0]
        y = mathieu_sem(m, q, v)[0]
        plt.plot(v, y, color=colors[i])
        leg.append(f"m = {m}")

    plt.title("Mathieu se2m+2")
    plt.xlabel("v")
    plt.ylabel("se")
    plt.legend(leg, loc="lower left")

    # ----------------------------------------------------------
    # Next: odd se functions (m = 1,3,5,7)
    plt.figure(2)
    leg = []

    for i, m in enumerate(range(1, 8, 2)):
        y = mathieu_sem(m, q, v)[0]
        plt.plot(v, y, color=colors[i])
        leg.append(f"m = {m}")

    plt.title("Mathieu se2m+1")
    plt.xlabel("v")
    plt.ylabel("se")
    plt.legend(leg, loc="lower left")

    # ----------------------------------------------------------
    # Plot se and derivative together
    N = 250
    u = np.linspace(0, 2 * np.pi, N)
    m = 3
    q = 5.0

    # Compute Mathieu function and derivative
    y, yd = mathieu_sem(m, q, u)

    txt = "m = %d\nq = %5.2f" % (m,q)
    plt.figure(3)
    plt.plot(u, y, label="se")
    plt.plot(u, yd, label="sed")
    plt.title("se and sed")
    plt.text(2.7,2.5,txt)
    plt.xlabel("u")
    plt.ylabel("se, sed")
    plt.legend()

    plt.show()

    
# ----------------------------------------------------------
# Run the function
if __name__ == "__main__":
    plot_mathieu_se()
