# ScipyMathieuTesting

This repo holds Python programs to plot and test the existing (and possibly
future) Mathieu functions in Scipy.  These tests were developed by
Stuart Brorson and students at Northeastern University in 2025 --
2026.  They highlight deficiencies in the existing Mathieu
implementations and serve as validation tests for the new Mathieu
implementation.  The figures in the paper under https://github.com/brorson/ScipyMathieuPaper/tree/main were created using stuff from this directory.

# Plotting functions:

* plot_eigs.py, plot_eigs_zoomed_out.py -- Plots Mathieu characteristic
values a, b.

* plot_mathieu_ce.py, plot_mathieu_se.py, plot_mathieu_modmc1.py,
plot_mathieu_modms1.py, plot_mathieu_modmc2.py, plot_mathieu_modms2.py
-- Plots Mathieu functions.

* plot_first_deriv_angfuncs.py, plot_first_deriv_angfuncs.py -- Plots
First derivative of Mathieu functions

# Selected test functions:

These were used to create the heat maps presented in the paper.

* simple_tests.py -- Creates two plots: For the existing Scipy Mathieu
impl the first plot shows the jumps in
the Mathieu eigenvalue for different order m and parameter q.  The
second plot shows how ce(m,q,z) jumps from one order to another for m
= 6 and q sweeps from 25.0 to 25.9.  This is a bug which is fixed in
the proposed new Mathieu impl.

* plot_round_trip_err_angfuncs.py -- Computes the second derivative
of the angular Matheu fcns for inputs m and q using a 6th order finite
difference formula.  Then subtracts the RHS of the Mathieu ODE to get
a residual.  Then scales the residual to get relative error.  Then
plots the rel error as a heat map vs. m and q.

* plot_round_trip_err_modfuncs.py -- Computes the second derivative
of the modified Matheu fcn for inputs m and q using a 6th order finite
difference formula.  Then subtracts the RHS of the Mathieu ODE to get
a residual.  Then scales the residual to get relative error.  Then
plots the rel error as a heat map vs. m and q.

* plot_wronskian_errs_modfcns.py -- Computes the Wronskian of pairs of
Mathieu fcns, then subtracts 2/pi which is the theoretical value.  The
result is a residual.  Then scales the residual to get relative error.
Then plots the rel error as a heat map vs. m and q.

* plot_wronskian_errs_bessels.py -- Similar to the above, computes the
  Wronskian of the J and Y Bessel functions.  This is just an accuracy
  check on the Bessel fcns provided by the Cephes library.

*  plot_first_deriv_err_angfuncs.py, plot_first_deriv_err_modfuncs.py --
Computes the first derivative of the fncs using a finite-difference
scheme, and then compares to the value returned by the impl.  Creates a
heat map of the error vs. m and q.

# Other files:

This directory holds many different numerical experiments in varying
states of completion.  Of import:

*  GVs -- Many files were used to perform golden-value testing.
Golden values were generated for the angular functions using a
collocation method, and then compared against the values returned by
the implementation.  This work is largely complete, but not organized.

