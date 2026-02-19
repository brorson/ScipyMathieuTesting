# SciPy Mathieu Function Testing

Validates SciPy's `mathieu_cem` and `mathieu_sem` against golden values (GVs) generated from MATLAB collocation, across a grid of orders `m` and parameter values `q`. This repository was created in Fall 2025.

## Repository Structure

```
.
├── code_to_generate_gvs/       # Scripts to generate golden value CSVs
│   ├── mp_math_ce.py           # Generate CE golden values using mpmath
│   ├── mp_math_se.py           # Generate SE golden values using mpmath
│   ├── np_gvs_ce.py            # Generate CE golden values using NumPy
│   └── np_gvs_se.py            # Generate SE golden values using NumPy
├── python_generated_gvs/       # Input golden value CSV files
├── python_not_passed_ce/       # Failure outputs for mathieu_cem
├── python_not_passed_se/       # Failure outputs for mathieu_sem
```

## Requirements

```bash
pip install numpy scipy pandas
```

## Usage

No path configuration needed — all paths are relative to the script location.

## Output

Failing test points are saved as CSVs with columns `m, v, gv, sp, diff`, where `gv` is the golden value, `sp` is SciPy's output, and `diff` is the absolute difference. Default tolerance: `atol = 2e-6`.