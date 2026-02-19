#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy.special import mathieu_cem, mathieu_sem
import glob
import os


# =======================
# CONFIG
# =======================

# Base directory: folder containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

BASE_DIR   = os.path.join(SCRIPT_DIR, "python_generated_gvs")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "python_not_passed")
ATOL = 2e-6


# =======================
# SETUP
# =======================

# Create output directory if missing
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Get ALL csv files
ALL_FILES = glob.glob(os.path.join(BASE_DIR, "*.csv"))

# Separate CE and SE
CE_FILES = [f for f in ALL_FILES if "mathieu_ce_gvs_q" in f]
SE_FILES = [f for f in ALL_FILES if "mathieu_se_gvs_q" in f]

print("Found", len(CE_FILES), "CE files")
print("Found", len(SE_FILES), "SE files")


# Row headers (v values)
row_headers = np.linspace(-np.pi, np.pi, 2500)
row_index = range(2500)


# =======================
# HELPER
# =======================

def extract_q(filepath):
    """
    Extract q value from filename.
    Example: mathieu_ce_gvs_q-30.000000.csv -> -30.0
    """
    name = os.path.basename(filepath)

    q_part = name.split("_q")[-1].replace(".csv", "")

    return float(q_part)


# =======================
# PROCESS CE FILES
# =======================

print("\nProcessing CE files...")

for file in sorted(CE_FILES):

    print("\nReading:", os.path.basename(file))

    data = np.loadtxt(file, delimiter=",", skiprows=1)

    q = extract_q(file)

    print("q =", q)

    diff_npass = []
    sp_npass = []
    gv_npass = []
    m_npass = []
    v_npass = []

    for v, v_ind in zip(row_headers, row_index):

        for m in range(1, 35):

            sp = mathieu_cem(m, q, v)[0]

            gv = data[v_ind, m + 1]   # CSV offset

            if not np.isclose(sp, gv, atol=ATOL):

                diff = abs(sp - gv)

                diff_npass.append(diff)
                sp_npass.append(sp)
                gv_npass.append(gv)
                m_npass.append(m)
                v_npass.append(v)

    results_df = pd.DataFrame({
        "m": m_npass,
        "v": v_npass,
        "gv": gv_npass,
        "sp": sp_npass,
        "diff": diff_npass
    })

    if not results_df.empty:

        out_file = f"{OUTPUT_DIR}/not_passed_ce_q{q}.csv"

        results_df.to_csv(out_file, index=False)

        print("Saved:", out_file)

    else:
        print("All CE tests passed for q =", q)


# =======================
# PROCESS SE FILES
# =======================

print("\nProcessing SE files...")

for file in sorted(SE_FILES):

    print("\nReading:", os.path.basename(file))

    data = np.loadtxt(file, delimiter=",", skiprows=1)

    q = extract_q(file)

    print("q =", q)

    diff_npass = []
    sp_npass = []
    gv_npass = []
    m_npass = []
    v_npass = []

    for v, v_ind in zip(row_headers, row_index):

        for m in range(1, 36):

            sp = mathieu_sem(m, q, v)[0]

            gv = data[v_ind, m]   # Original indexing

            if not np.isclose(sp, gv, atol=ATOL):

                diff = abs(sp - gv)

                diff_npass.append(diff)
                sp_npass.append(sp)
                gv_npass.append(gv)
                m_npass.append(m)
                v_npass.append(v)

    results_df = pd.DataFrame({
        "m": m_npass,
        "v": v_npass,
        "gv": gv_npass,
        "sp": sp_npass,
        "diff": diff_npass
    })

    if not results_df.empty:

        out_file = f"{OUTPUT_DIR}/not_passed_se_q{q}.csv"

        results_df.to_csv(out_file, index=False)

        print("Saved:", out_file)

    else:
        print("All SE tests passed for q =", q)


# =======================
# DONE
# =======================

print("\nFinished running all tests.")

# Quick sanity check
print("Sanity check:", mathieu_sem(1, -0.001, -180))
