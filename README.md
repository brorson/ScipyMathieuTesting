# ScipyMathieuTesting

This repo holds Python programs to test the existing (and possibly
future) Mathieu functions in Scipy.  These tests were developed by
Stuart Brorson and students at Northeastern University in 2025 --
2026.  They highlight deficiencies in the existing Mathieu
implementations and serve as validation tests for the new Mathieu
implementation. 

# The following tests are of interest:

* simple_tests.py -- Creates two plots: For the existing Scipy Mathieu
impl the first plot shows the jumps in
the Mathieu eigenvalue for different order m and parameter q.  The
second plot shows how ce(m,q,z) jumps from one order to another for m
= 6 and q sweeps from 25.0 to 25.9.  This is a bug which is fixed in
the proposed new Mathieu impl.

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



