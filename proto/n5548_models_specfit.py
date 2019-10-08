# Model for n5548 in new format
#
# This was further modified so it can be directly ingested by
# the specutils' fit_lines function. This function does not
# seem to accept compound models made out of subclasses of the
# original astropy.modelling.models classes, as defined in
# n5548_models.py. It has to be build out of classes taken
# directly from astropy.modeling.models, and not subclassed
# from Fittable1DModel.

#
# This model differs from the original model specified in
# n5548_models.py, by lacking a ccmext extinction component.
# It also lacks the bounds and ties associated with each
# component's parameters.

from astropy import units as u
from astropy.units import Quantity
import astropy.modeling.models as models


# convert km/s to wavelength range in Angstrom

def kms2Angstrom(value, x0):
    return (Quantity(value *u.km/u.s).to(u.Angstrom, equivalencies=u.doppler_optical(x0*u.Angstrom))).value - x0


# convert fwhm in km/s (as in Kriss's model specs)
# to stddedv in Angstrom (as in specutils/astropy)

def fwhm2stddev(fwhm, mean):
    return kms2Angstrom(fwhm, mean) / 2.355


# convert flux (as in Kriss's model specs)
# to amplitude (as in specutils/astropy)

def flux2amplitude(flux, mean, fwhm):
    return flux / kms2Angstrom(fwhm, mean) * 0.937415


model1 = \
    models.PowerLaw1D(
             amplitude = 6.586200E-14,
             x_0 = 1000.0,
             alpha = 0.4819233
             ) \
+ \
    models.Gaussian1D(
             amplitude = flux2amplitude(2.000000E-14, 1195.006, 861.4926),
             mean = 1195.006,
             stddev = fwhm2stddev(861.4926, 1195.006)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(1.438015E-14, 1226.392, 861.4926),
             mean = 1226.392,
             stddev = fwhm2stddev(861.4926, 1226.392)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(2.020000E-13, 1236.729, 255.4998),
             mean = 1236.729,
             stddev = fwhm2stddev(255.4998, 1236.729)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(5.474183E-13, 1235.996, 861.4926),
             mean = 1235.996,
             stddev = fwhm2stddev(861.4926, 1235.996)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(3.948799E-12, 1235.138, 3040.59),
             mean = 1235.138,
             stddev = fwhm2stddev(3040.59, 1235.138)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(2.156964E-12, 1231.589, 8133.099),
             mean = 1231.589,
             stddev = fwhm2stddev(8133.099, 1231.589)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(0.536853E-12, 1237.643, 18183.71),
             mean = 1237.643,
             stddev = fwhm2stddev(18183.71, 1237.643)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(1.217935E-14, 1259.753, 255.4998),
             mean = 1259.753,
             stddev = fwhm2stddev(255.4998, 1259.753)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(1.217935E-14, 1263.803, 255.4998),
             mean = 1263.803,
             stddev = fwhm2stddev(255.4998, 1263.803)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(6.219548E-15, 1259.533, 861.4926),
             mean = 1259.533,
             stddev = fwhm2stddev(861.4926, 1259.533)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(6.219548E-15, 1263.582, 861.4926),
             mean = 1263.582,
             stddev = fwhm2stddev(861.4926, 1263.582)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(0.221692E-13, 1258.659, 3040.59),
             mean = 1258.659,
             stddev = fwhm2stddev(3040.59, 1258.659)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(2.221692E-13, 1262.705, 3040.59),
             mean = 1262.705,
             stddev = fwhm2stddev(3040.59, 1262.705)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(3.185217E-13, 1255.042, 8133.099),
             mean = 1255.042,
             stddev = fwhm2stddev(8133.099, 1255.042)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(1.185217E-13, 1259.077, 8133.099),
             mean = 1259.077,
             stddev = fwhm2stddev(8133.099, 1259.077)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(2.287816E-13, 1263.24, 18183.71),
             mean = 1263.24,
             stddev = fwhm2stddev(18183.71, 1263.24)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(-1.000000E-15, 1194., 3683.102),
             mean = 1194.,
             stddev = fwhm2stddev(3683.102, 1194.)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(-5.349239E-13, 1235., 3683.102),
             mean = 1235.,
             stddev = fwhm2stddev(3683.102, 1235.)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(-10.921041E-14, 1258., 3683.102),
             mean = 1258.,
             stddev = fwhm2stddev(3683.102, 1258.)
             ) \
+ \
    models.Gaussian1D(
             amplitude=flux2amplitude(-8.921041E-14, 1262.044, 3683.102),
             mean = 1262.044,
             stddev = fwhm2stddev(3683.102, 1262.044)
             )