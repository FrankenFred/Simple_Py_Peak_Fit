# Simple_Py_Peak_Fit
Simple Python 3x script that reads in a spectrum and uses HITRAN line data to perform a Voight profile NLLSQ fit (with offset). Change line, experimental and file properties as needed. Available in Jupyter notebook and .py format.

## Version 1.0.0-beta:
Current version estimates peaks before fitting and compares to the sample spectrum. This is to give the user the ability to optimize the initial guesses for the fitting routine. If the peak centres are off by too much, then the fit does not converge properly.
* Reordered code into sections
* Fits four lines in mid IR in Abs vs current (mA)

## Version 0:
Preliminary version shows fitting capability and is hard coded to:
* Fit a single line in the MIR (used experimentally)
* Only fits in the Absorption/wavenumber domain

