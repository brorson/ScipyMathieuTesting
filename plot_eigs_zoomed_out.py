import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_a, mathieu_b


def plot_eigs_zoomed():
    """
    Make plot of Mathieu eigenvalues vs. q (zoomed out view)
    Extended range: q from -100 to 100, eigenvalues from -100 to 300
    """
    # Number of sample points
    N = 1000

    # Extended domain of q values (including negative q)
    qs = np.linspace(-100, 100, N)

    # More eigenvalues to track for the zoomed out view
    Ne = 20

    # Calculate a eigenvalues (ce eigs)
    a_vals = np.zeros((len(qs), Ne))
    print('Calculating a eigenvalues ...')

    for m in range(Ne):
        for i, q in enumerate(qs):
            a_vals[i, m] = mathieu_a(m, q)

    # Make plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot a eigenvalues (ce eigs)
    for j in range(Ne):
        if j == 0:
            h1, = ax.plot(qs, a_vals[:, j], 'b-', linewidth=1.5, label='ce eigs')
        else:
            ax.plot(qs, a_vals[:, j], 'b-', linewidth=1.5)

    # Calculate b eigenvalues (se eigs)
    b_vals = np.zeros((len(qs), Ne))
    print('Calculating b eigenvalues ...')

    for m in range(1, Ne + 1):
        for i, q in enumerate(qs):
            b_vals[i, m - 1] = mathieu_b(m, q)

    # Plot b eigenvalues (se eigs)
    for j in range(Ne):
        if j == 0:
            h2, = ax.plot(qs, b_vals[:, j], 'r--', linewidth=1.5, label='se eigs')
        else:
            ax.plot(qs, b_vals[:, j], 'r--', linewidth=1.5)

    # Set axis limits for zoomed out view
    ax.set_xlim([-100, 100])
    ax.set_ylim([-100, 300])

    ax.set_title('First Mathieu eigenvalues vs. q')
    ax.set_xlabel('q')
    ax.set_ylabel('eigenvalue')
    ax.legend(loc='upper left')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_eigs_zoomed()