# Simple_Py_Peak_Fit
Simple Python 3x script that reads in a spectrum and uses HITRAN line data to perform a Voight profile NLLSQ fit (with offset). Change line, experimental and file properties as needed. Available in Jupyter notebook and .py format.

## Version 1.1.0:
Added functionality to fit raw spectrum in Volts vs Current. Fits a polynomial as the BKG spectrum (I0). Polynomial is procedural and so the order can be specified by handing the function more arguments. Data needs to be trimmed as fit function will fit the entire visible spectrum in the peak estimation phase.

## Version 1.0.0-beta:
Current version estimates peaks before fitting and compares to the sample spectrum. This is to give the user the ability to optimize the initial guesses for the fitting routine. If the peak centres are off by too much, then the fit does not converge properly.
* Reordered code into sections
* Fits four lines in mid IR in Abs vs current (mA)

## Version 0:
Preliminary version shows fitting capability and is hard coded to:
* Fit a single line in the MIR (used experimentally)
* Only fits in the Absorption/wavenumber domain

