# Model for n5548 in new format
#
# This was further modified so it can be directly ingested by
# the specutils' fit_lines function. This function does not
# accept compound models made out of subclasses of the
# original astropy.modelling.models classes, as defined in
# n5548_models.py

from astropy import units as u
from astropy.units import Quantity
import astropy.modeling.models as models


# TODO norm must be converted to amplitude


# function for converting stddev from km/s (as in Kriss's model
# specs) to Angstrom (as in specutils/astropy)

def stddev_function(stddev, mean):
    return (Quantity(stddev * u.km / u.s).to(u.Angstrom, equivalencies=u.doppler_optical(mean * u.Angstrom)) - (mean * u.Angstrom)).value


model1 = \
    models.PowerLaw1D(
             amplitude = 6.586200E-14,
             x_0 = 1000.0,
             alpha = 0.4819233
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.000000E-14,
             mean = 1195.006,
             stddev = stddev_function(861.4926, 1195.006)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.438015E-14,
             mean = 1226.392,
             stddev = stddev_function(861.4926, 1226.392)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.020000E-13,
             mean = 1236.729,
             stddev = stddev_function(255.4998, 1236.729)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 5.474183E-13,
             mean = 1235.996,
             stddev = stddev_function(861.4926, 1235.996)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 3.948799E-12,
             mean = 1235.138,
             stddev = stddev_function(3040.59, 1235.138)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.156964E-12,
             mean = 1231.589,
             stddev = stddev_function(8133.099, 1231.589)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 0.536853E-12,
             mean = 1237.643,
             stddev = stddev_function(18183.71, 1237.643)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.217935E-14,
             mean = 1259.753,
             stddev = stddev_function(255.4998, 1259.753)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.217935E-14,
             mean = 1263.803,
             stddev = stddev_function(255.4998, 1263.803)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 6.219548E-15,
             mean = 1259.533,
             stddev = stddev_function(861.4926, 1259.533)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 6.219548E-15,
             mean = 1263.582,
             stddev = stddev_function(861.4926, 1263.582)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 0.221692E-13,
             mean = 1258.659,
             stddev = stddev_function(3040.59, 1258.659)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.221692E-13,
             mean = 1262.705,
             stddev = stddev_function(3040.59, 1262.705)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 3.185217E-13,
             mean = 1255.042,
             stddev = stddev_function(8133.099, 1255.042)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.185217E-13,
             mean = 1259.077,
             stddev = stddev_function(8133.099, 1259.077)
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.287816E-13,
             mean = 1263.24,
             stddev = stddev_function(18183.71, 1263.24)
             ) \
+ \
    models.Gaussian1D(
             amplitude = -1.000000E-15,
             mean = 1194.,
             stddev = stddev_function(3683.102, 1194.)
             ) \
+ \
    models.Gaussian1D(
             amplitude = -5.349239E-13,
             mean = 1235.,
             stddev = stddev_function(3683.102, 1235.)
             ) \
+ \
    models.Gaussian1D(
             amplitude = -10.921041E-14,
             mean = 1258.,
             stddev = stddev_function(3683.102, 1258.)
             ) \
+ \
    models.Gaussian1D(
             amplitude = -8.921041E-14,
             mean = 1262.044,
             stddev = stddev_function(3683.102, 1262.044)
             )