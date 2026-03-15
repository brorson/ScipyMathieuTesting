import numpy as np
import scipy as sp
from scipy.special import *

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# ============================================================
# This script creates plots of the error of the modified
# Mathieu functions using Wronskian deviation.
# ============================================================

# Parameters
ms = np.arange(1, 51)
qs = np.logspace(-4, 4, 30)

# Domain
N = 100
v = np.linspace(2, 10, N)[:, None]  # Column vector

# True value of Wronskian
wtrue = 2 / np.pi

# ============================================================
# Heatmap configuration (mirrors test_gvs_ce.py)
# ============================================================

# Colors (same palette as test_gvs_ce.py)
PASSED_COLOR     = "#8BCF8B"   # light green
NOT_PASSED_COLOR = "#C44E52"   # medium red
OTHER_COLOR      = "#F2F2F2"   # neutral grey

# A point is considered "passed" when its log10(rms error) is below this.
PASS_THRESHOLD = 2e-4           # i.e. rms error < 1e-6


# ============================================================
# Helper: contourf plot (original)
# ============================================================

def make_plot(fig_num, X, Y, errs, title):
    plt.figure(fig_num, figsize=(8, 6))
    levels = np.arange(-20, 5, 5)
    cs = plt.contourf(X, Y, errs, levels=levels, extend='both')
    plt.clabel(cs, inline=True)
    plt.xlabel("Order m")
    plt.ylabel("log10(q)")
    plt.title(title)
    plt.clim(-20, 10)
    plt.colorbar()
    plt.show()


# ============================================================
# Helper: continuous error heatmap
# ============================================================

def save_error_heatmap(errs, m_values, log10_q_values, title):
    fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True)

    heat = errs.T                  # shape: (n_q, n_m)

    extent = [
        m_values[0] - 0.5, m_values[-1] + 0.5,
        log10_q_values[0],  log10_q_values[-1],
    ]

    im = ax.imshow(
        heat,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="RdYlGn_r",
        vmin=-20,
        vmax=0,
        extent=extent,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10(rms Wronskian error)")

    ax.set_xlabel("Order m")
    ax.set_ylabel("log10(q)")
    ax.set_title(title)
    ax.set_box_aspect(1)

    xtick_step = max(1, len(m_values) // 12)
    ax.set_xticks(m_values[::xtick_step])
    ax.set_ylim(log10_q_values[0], log10_q_values[-1])

    plt.show()


# ============================================================
# Helper: binary pass/fail heatmap
# ============================================================

def save_binary_heatmap(mask, m_values, log10_q_values,
                        active_label, active_color, title):
    fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True)

    heat = mask.astype(int).T      # shape: (n_q, n_m)
    cmap = ListedColormap([OTHER_COLOR, active_color])
    extent = [
        m_values[0] - 0.5, m_values[-1] + 0.5,
        log10_q_values[0],  log10_q_values[-1],
    ]

    im = ax.imshow(
        heat,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1,
        extent=extent,
    )

    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["other", active_label])

    ax.set_xlabel("Order m")
    ax.set_ylabel("log10(q)")
    ax.set_title(title)
    ax.set_box_aspect(1)

    xtick_step = max(1, len(m_values) // 12)
    ax.set_xticks(m_values[::xtick_step])
    ax.set_ylim(log10_q_values[0], log10_q_values[-1])

    plt.show()


# ============================================================
# Helper: combined pass + not-passed heatmap
# ============================================================

def save_combined_heatmap(pass_mask, m_values, log10_q_values, title):
    fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True)

    heat = pass_mask.astype(int).T  # 0 = not-passed, 1 = passed
    cmap = ListedColormap([NOT_PASSED_COLOR, PASSED_COLOR])
    extent = [
        m_values[0] - 0.5, m_values[-1] + 0.5,
        log10_q_values[0],  log10_q_values[-1],
    ]

    im = ax.imshow(
        heat,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1,
        extent=extent,
    )

    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["not-passed", "passed"])

    ax.set_xlabel("Order m")
    ax.set_ylabel("log10(q)")
    ax.set_title(title)
    ax.set_box_aspect(1)

    xtick_step = max(1, len(m_values) // 12)
    ax.set_xticks(m_values[::xtick_step])
    ax.set_ylim(log10_q_values[0], log10_q_values[-1])

    plt.show()


# ============================================================
# Shared helper: generate all three heatmaps for a pair of functions
# ============================================================

def generate_heatmaps(errs, label):
    log10_q = np.log10(qs)
    pass_mask = errs < PASS_THRESHOLD   # shape: (n_m, n_q)

    n_pass     = int(pass_mask.sum())
    n_total    = pass_mask.size
    n_not_pass = n_total - n_pass

    # 1. Continuous error heatmap (always produced)
    save_error_heatmap(
        errs, ms, log10_q,
        title=f"log10(rms Wronskian error) — {label}",
    )

    # 2–4. Pass/fail heatmaps (only when both classes are present)
    if n_pass == 0 or n_not_pass == 0:
        only = "passed" if n_pass == n_total else "not-passed"
        print(f"  Pass/fail heatmaps skipped for {label} (all points are {only}).")
        return

    save_binary_heatmap(
        pass_mask, ms, log10_q,
        active_label="passed",
        active_color=PASSED_COLOR,
        title=f"Passed heatmap — {label} (threshold: log10 err < {PASS_THRESHOLD})",
    )
    save_binary_heatmap(
        ~pass_mask, ms, log10_q,
        active_label="not-passed",
        active_color=NOT_PASSED_COLOR,
        title=f"Not-passed heatmap — {label} (threshold: log10 err < {PASS_THRESHOLD})",
    )
    save_combined_heatmap(
        pass_mask, ms, log10_q,
        title=f"Passed + not-passed heatmap — {label}",
    )


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
        rms = np.sqrt(np.mean((w - wtrue)**2))
        errs[i, j] = np.log10(rms)
        X[i, j] = m
        Y[i, j] = np.log10(q)

make_plot(1, X, Y, errs,
          "Log10 of Wronskian rms error -- Mc1 Mc2")

generate_heatmaps(errs, "Mc1_Mc2")


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
        rms = np.sqrt(np.mean((w - wtrue)**2))
        errs[i, j] = np.log10(rms)
        X[i, j] = m
        Y[i, j] = np.log10(q)

make_plot(2, X, Y, errs,
          "Log10 of Wronskian rms error -- Ms1 Ms2")

generate_heatmaps(errs, "Ms1_Ms2")


# ============================================================
# Mc1 & Ms2  — excluded (different eigenvalues, invalid Wronskian test)
# ============================================================
# See SDB comment in original script.


# ============================================================
# Ms1 & Mc2  — excluded (different eigenvalues, invalid Wronskian test)
# ============================================================
# See SDB comment in original script.


print("\nAll plots generated successfully.")
input("Press Enter to close...")