# SciPy Mathieu Function Testing - Golden Values

Validates SciPy's `mathieu_cem` and `mathieu_sem` against golden values (GVs) generated from MATLAB collocation, across a grid of orders `m` and parameter values `q`. This repository was created in Spring 2026.

## Repository Structure

```
.
├── Code_to_generate_gvs/   # Scripts to generate golden value CSVs
├── python_generated_gvs/   # Input golden value CSV files
├── python_not_passed/      # Combined CE + SE failure outputs
├── python_not_passed_ce/   # Failure outputs for mathieu_cem
├── python_not_passed_se/   # Failure outputs for mathieu_sem
├── mathieu_compare_gvs.py  # Combined CE + SE comparison script
├── test_gvs_ce.py          # Standalone test for mathieu_cem
└── test_gvs_se.py          # Standalone test for mathieu_sem
```

## Requirements

```bash
pip install numpy scipy pandas
```

## Usage

```bash
python test_gvs.py  # Run both CE and SE tests
python test_gvs_ce.py          # CE only
python test_gvs_se.py          # SE only
```

No path configuration needed — all paths are relative to the script location.

## Output

Failing test points are saved as CSVs with columns `m, v, gv, sp, diff`, where `gv` is the golden value, `sp` is SciPy's output, and `diff` is the absolute difference. Default tolerance: `atol = 2e-6`.