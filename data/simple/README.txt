Version 1 -- 2015-Apr-06 -- Jeff Valenti

This directory contains simple test cases that can be used to assess
the numerical behavior of spectrum fitting software. For noiseless
cases, IDL recovers the input parameters to high-precision. For
noisy cases, other solvers may legitimately disagree on the optimal
solution.

The directory contains three IDL programs:

func.pro -- given an input wavelength vector and parameter vector,
calculate a noiseless synthetic spectrum that is the sum of one or
more Gaussians with no continuum.

obs.pro -- calculate four simulated spectra, two with one Gaussian
and two with two Gaussians. For each case, one simulated spectrum
has no noise and the other has noise.

fit.pro -- fit each of the simulated spectra using the MPFITFUN
routine written by Craig Markwardt. For each of the four simulated
spectra try three separate fits with different fixed parameters.
First fix Gaussian center(s) and sigma(s). Next only fix sigma(s).
Finally, don't fix any parameters. Write a log file and produce
postscript plots.

Only the IDL routines mpfit.pro and mpfitfun.pro are required to
run this software. There are no other IDL dependencies. To run the
software, untar the directory, start IDL, and the run 'fit' at the
command prompt. All output files (except the PNG versions of the
postscript) will be regenerated.

The four simulated observations are written to four text files
with a header (delimited by #) and then two or three columns,
depending on whether the observations are noiseless or include
noise. For noiseless spectra, column 1 has wavelength (in microns)
and column 2 has flux (in arbitrary units). Noisy spectra have
a third column, which is the flux uncertainy (in the same arbitrary
units as flux). The random number generator seed is specified, so
rerunning the code yields identical noise. The four observation
files are:

obs1g.dat -- one Gaussian, no noise
obs1gn.dat -- one Gaussian with noise
obs2g.dat -- two overlapping Gaussians, no noise
obs2gn.dat -- two overlapping Gaussian with noise

Running the IDL program 'fit' (re)generates the one log file and
two postscript files:

fit.log -- summary of all 12 fitting examples (four input spectra
times three flavors of fixed parameters).

obs1gn.ps -- plot of the fitting example with one noisy Gaussian
and no fixed parameters.

obs2gn.ps -- plot of the fitting example with two noisy Gaussians
and no fixed parameters.

The log file has the following columns:

Case -- indicates model parameters or name of observation fitted
It -- Number of iternations to converge (0 for model description)
S -- Exit status from MPFITFUNC (see documentation if you care)
Fixed -- Flags indicating which of the parameters were fixed (1)
Row -- Type of data in the row as follows:
         Mod -- assumed model parameters
         Ini -- initial parameter guesses fed to solver
         Fin -- final parameter values returned by solver
         Unc -- parameter uncertainties returned by solver
         Cor -- Parameter correlation matrix returned by solver
Amp1 -- amplitude of Gaussian 1 [arbitrary units]
Cen1 -- central wavelength of Gaussian 1 [in microns]
Sig1 -- standard deviation of Gaussian 1 [in microns]
Amp2 -- amplitude of Gaussian 2 [arbitrary units]
Cen2 -- central wavelength of Gaussian 2 [in microns]
Sig2 -- standard deviation of Gaussian 2 [in microns]

The parameter correlation matrix Cor = Cov / (Unc # Unc), where #
denotes an outer product. It is good to renormalize the coveriance
matrix in this way because the uncertainty in Amp is three orders
of magnitude larger than the uncertainty in Cen and Sig.

The plots are fairly straightforward. The black X symbols are the
noisy observation. The blue curve is the noiseless model. The red
curve is the model fitted to the noisy observations. For the two
Gaussian case, the two model Gaussians are also shown separately
in green.
