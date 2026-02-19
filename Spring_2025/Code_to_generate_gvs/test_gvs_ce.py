#!/usr/bin/env python
"""
test_gvs_ce.py
==============
Test SciPy's mathieu_cem against golden values (GVs).

Golden values are generated externally (e.g. from MATLAB collocation)
and stored as CSV files. This script reads each GV file, compares
SciPy's mathieu_cem output at each (m, q, v) grid point, and writes
out any non-passing results.

Directory layout expected:
    <GV_INPUT_DIR>/mathieu_ce_gvs_q<qval>.csv   -- input GV files
    <RESULTS_OUTPUT_DIR>/not_passed_ce_q<qval>.csv -- output failures

CSV format of GV files:
    - First row: header (skipped)
    - Columns: v, ce(m=0), ce(m=1), ..., ce(m=M_MAX)
    - Rows: 2500 points from v = -pi to pi

Usage:
    python test_gvs_ce.py

Adjust the configuration constants below as needed.
"""

import numpy as np
import pandas as pd
import glob
import os
import sys

from scipy.special import mathieu_cem

# ============================================================
# Configuration -- adjust these paths and parameters as needed
# ============================================================
# '/Users/kushala/ScipyMathieuTesting_testing/python_generated_gvs/mathieu_ce_gvs_q-1.000000.csv'
# Directory containing golden value CSV files
# Base directory: folder containing this script (portable across machines)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Directory containing golden value CSV files
GV_INPUT_DIR = os.path.join(SCRIPT_DIR, "python_generated_gvs")

# Glob pattern to match ce golden value files
GV_FILE_PATTERN = os.path.join(GV_INPUT_DIR, "mathieu_se_gvs_q-0.001000*")

# Directory to write non-passing results
RESULTS_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "python_not_passed_ce")

# Number of evaluation points in v (from -pi to pi)
N_POINTS = 2500

# Range of orders m to test: m = M_START, M_START+1, ..., M_END-1
M_START = 1
M_END = 35  # exclusive upper bound, so tests m=1..34

# Column offset: in the GV CSV, ce(m) is at column index (m + COL_OFFSET)
# This accounts for the first column being the v values
COL_OFFSET = 1

# Absolute tolerance for np.isclose comparison
ATOL = 2e-6

# ============================================================
# Helper functions
# ============================================================

def extract_q_from_filename(filepath):
    """
    Extract the q parameter value from the GV filename.

    Expects filenames like: .../mathieu_ce_gvs_q1.0.csv
    Finds the substring after the last 'q' before '.csv'.
    """
    basename = os.path.basename(filepath)
    # Find 'q' and extract the number between 'q' and '.csv'
    q_start = basename.rfind('q') + 1
    q_str = basename[q_start:].replace('.csv', '')
    return float(q_str)


def load_gv_data(filepath):
    """Load golden value CSV, skipping the header row."""
    return np.loadtxt(filepath, delimiter=',', skiprows=1)


def test_ce_for_file(filepath, v_grid):
    """
    Compare SciPy mathieu_cem against golden values for a single q.

    Parameters
    ----------
    filepath : str
        Path to the golden value CSV file.
    v_grid : np.ndarray
        Array of v values (angles) at which to evaluate.

    Returns
    -------
    results_df : pd.DataFrame
        DataFrame of non-passing test points with columns:
        m, v, gv, sp, diff
    q : float
        The q parameter extracted from the filename.
    """
    q = extract_q_from_filename(filepath)
    data = load_gv_data(filepath)

    diff_npass = []
    sp_npass = []
    gv_npass = []
    m_npass = []
    v_npass = []

    for v_ind, v in enumerate(v_grid):
        for m in range(M_START, M_END):
            sp_val = mathieu_cem(m, q, v)[0]
            gv_val = data[v_ind, m + COL_OFFSET]

            if not np.isclose(sp_val, gv_val, atol=ATOL):
                diff = np.abs(sp_val - gv_val)
                diff_npass.append(diff)
                sp_npass.append(sp_val)
                gv_npass.append(gv_val)
                m_npass.append(m)
                v_npass.append(v)

    results_df = pd.DataFrame({
        'm': m_npass,
        'v': v_npass,
        'gv': gv_npass,
        'sp': sp_npass,
        'diff': diff_npass
    })

    return results_df, q


# ============================================================
# Main
# ============================================================

def main():
    # Create output directory if it doesn't exist
    os.makedirs(RESULTS_OUTPUT_DIR, exist_ok=True)

    # Build the v grid: N_POINTS equally spaced from -pi to pi
    v_grid = np.linspace(-np.pi, np.pi, N_POINTS)

    # Find all matching GV files
    matching_files = sorted(glob.glob(GV_FILE_PATTERN))

    if not matching_files:
        print(f"No golden value files found matching: {GV_FILE_PATTERN}")
        print("Please check that the GV CSV files are in the correct directory.")
        sys.exit(1)

    print(f"Found {len(matching_files)} golden value file(s) for ce.")
    print(f"Testing orders m = {M_START} to {M_END - 1}")
    print(f"Tolerance: atol = {ATOL}")
    print(f"Grid: {N_POINTS} points from -pi to pi")
    print("-" * 60)

    total_failures = 0

    for filepath in matching_files:
        print(f"Processing: {filepath}")

        try:
            results_df, q = test_ce_for_file(filepath, v_grid)
        except Exception as e:
            print(f"  ERROR processing file: {e}")
            continue

        n_failures = len(results_df)
        total_failures += n_failures

        # Write results
        output_file = os.path.join(
            RESULTS_OUTPUT_DIR, f"not_passed_ce_q{q}.csv"
        )
        results_df.to_csv(output_file, index=False)

        if n_failures == 0:
            print(f"  q = {q}: ALL PASSED")
        else:
            print(f"  q = {q}: {n_failures} failures written to {output_file}")
            # Print summary statistics for failures
            if n_failures > 0:
                print(f"    Max diff: {results_df['diff'].max():.6e}")
                print(f"    Mean diff: {results_df['diff'].mean():.6e}")
                print(f"    Orders with failures: "
                      f"{sorted(results_df['m'].unique().tolist())}")

    print("-" * 60)
    print(f"Finished. Total failures across all q values: {total_failures}")


if __name__ == "__main__":
    main()
