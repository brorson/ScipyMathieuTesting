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

    plt.title("Mathieu se2n")
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

    plt.title("Mathieu se2n+1")
    plt.xlabel("v")
    plt.ylabel("se")
    plt.legend(leg, loc="lower left")

    plt.show()


# ----------------------------------------------------------
# Run the function
if __name__ == "__main__":
    plot_mathieu_se()