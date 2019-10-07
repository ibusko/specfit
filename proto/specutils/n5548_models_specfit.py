# Model for n5548 in new format
#
# This was further modified so it can be directly ingested by
# the specutils' fit_lines function. This function does not
# accept compound models made out of subclasses of the
# original astropy.modelling.models classes, as defined in
# n5548_models.py

import astropy.units as u
from astropy.units import Unit, Quantity
import astropy.modeling.models as models
import astropy.units.equivalencies as eq

# template for converting stddev from km/s (as in Kriss's model
# specs) to Angstrom (as in specutils/astropy)
#
# (STDDEV *u.km/u.s).to(u.Angstrom,equivalencies=u.doppler_optical(MEAN *u.Angstrom))-(MEAN *u.Angstrom)

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
             stddev = Quantity(861.4926*Unit('km/s')).to(Unit('Angstrom'),equivalencies='spectral')
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.438015E-14,
             mean = 1226.392,
             stddev = 861.4926
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.020000E-13,
             mean = 1236.729,
             stddev = 255.4998
             ) \
+ \
    models.Gaussian1D(
             amplitude = 5.474183E-13,
             mean = 1235.996,
             stddev = 861.4926
             ) \
+ \
    models.Gaussian1D(
             amplitude = 3.948799E-12,
             mean = 1235.138,
             stddev = 3040.59
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.156964E-12,
             mean = 1231.589,
             stddev = 8133.099
             ) \
+ \
    models.Gaussian1D(
             amplitude = 0.536853E-12,
             mean = 1237.643,
             stddev = 18183.71
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.217935E-14,
             mean = 1259.753,
             stddev = 255.4998
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.217935E-14,
             mean = 1263.803,
             stddev = 255.4998
             ) \
+ \
    models.Gaussian1D(
             amplitude = 6.219548E-15,
             mean = 1259.533,
             stddev = 861.4926
             ) \
+ \
    models.Gaussian1D(
             amplitude = 6.219548E-15,
             mean = 1263.582,
             stddev = 861.4926
             ) \
+ \
    models.Gaussian1D(
             amplitude = 0.221692E-13,
             mean = 1258.659,
             stddev = 3040.59
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.221692E-13,
             mean = 1262.705,
             stddev = 3040.59
             ) \
+ \
    models.Gaussian1D(
             amplitude = 3.185217E-13,
             mean = 1255.042,
             stddev = 8133.099
             ) \
+ \
    models.Gaussian1D(
             amplitude = 1.185217E-13,
             mean = 1259.077,
             stddev = 8133.099
             ) \
+ \
    models.Gaussian1D(
             amplitude = 2.287816E-13,
             mean = 1263.24,
             stddev = 18183.71
             ) \
+ \
    models.Gaussian1D(
             amplitude = -1.000000E-15,
             mean = 1194.,
             stddev = 3683.102
             ) \
+ \
    models.Gaussian1D(
             amplitude = -5.349239E-13,
             mean = 1235.,
             stddev = 3683.102
             ) \
+ \
    models.Gaussian1D(
             amplitude = -10.921041E-14,
             mean = 1258.,
             stddev = 3683.102
             ) \
+ \
    models.Gaussian1D(
             amplitude = -8.921041E-14,
             mean = 1262.044,
             stddev = 3683.102
             )