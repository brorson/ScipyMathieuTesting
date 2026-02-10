import numpy as np
import matplotlib.pyplot as plt
from scipy.special import mathieu_a, mathieu_b


def plot_eigs():
    """
    Make plot of Mathieu eigenvalues vs. q to reproduce plot
    on https://dlmf.nist.gov/28.2
    """
    # Number of sample points
    N = 1000

    # Domain of q values to examine (for plotting)
    qs = np.linspace(0, 10, N)

    # Number of each type of eigenvalue to track
    Ne = 6  # Ne of a and Ne of b

    # First plot ce eigs
    # Preallocate array to store values
    a_vals = np.zeros((len(qs), Ne))
    print('Calculating a eigenvalues ...')

    # Loop over order m
    for m in range(Ne):
        # Loop over qs
        for i, q in enumerate(qs):
            a_vals[i, m] = mathieu_a(m, q)

    # Make plot of eigenvalues vs. q
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot a eigenvalues (ce eigs)
    for j in range(Ne):
        if j == 0:
            h1, = ax.plot(qs, a_vals[:, j], 'b-', linewidth=2, label='ce eigs')
        else:
            ax.plot(qs, a_vals[:, j], 'b-', linewidth=2)

    # Next plot se eigs
    # Preallocate b array to store values
    b_vals = np.zeros((len(qs), Ne))

    print('Calculating b eigenvalues ...')

    # Loop over order m (b starts at m=1)
    for m in range(1, Ne + 1):
        # Loop over qs
        for i, q in enumerate(qs):
            b_vals[i, m - 1] = mathieu_b(m, q)

    # Plot b eigenvalues (se eigs)
    for j in range(Ne):
        if j == 0:
            h2, = ax.plot(qs, b_vals[:, j], 'r--', linewidth=2, label='se eigs')
        else:
            ax.plot(qs, b_vals[:, j], 'r--', linewidth=2)

    # Set axis limits to reproduce the DLMF plot
    ax.set_xlim([0, 10])
    ax.set_ylim([-5, 20])

    ax.set_title('First Mathieu eigenvalues vs. q')
    ax.set_xlabel('q')
    ax.set_ylabel('eigenvalue')
    ax.legend(loc='upper left')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_eigs()
