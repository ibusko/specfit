# Model for n5548 in new format

from custom_models import gaussian, powerlaw, ccmext

model1 = \
    powerlaw(name = 'powerlaw1',
             amp =   6.586200E-14,
             x_0 =   1000.0,
             alpha = 0.4819233,
             bounds = {'amp':   (0., 1.00E-11),
                       'x_0':   (0., 1.00E-11),
                       'alpha': (-5., 5.)},
             fixed = {'x_0': True}
             ) \
* \
    ccmext(name = 'extinction',
           ebmv = 0.01713,
           rv = 3.1,
           bounds = {'ebmv': (0., 1.),
                     'rv':   (3.1, 3.1)},
           fixed = {'ebmv': True,
                    'rv': True}
             ) \
+ \
    gaussian(name = 'C III 1176',
             norm = 2.000000E-14,
             mean = 1195.006,
             fwhm = 861.4926,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (1000., 2000.),
                       'skew': (1., 1.)},
             fixed = {'norm': True,
                      'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'Si III 1206',
             norm = 1.438015E-14,
             mean = 1226.392,
             fwhm = 861.4926,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (500., 2000.),
                       'skew': (1., 1.)},
             tied = {'fwhm': lambda m: 1.0 * m['C III 1176'].fwhm},
             fixed = {'norm': True,
                      'mean': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'Ly alpha - NLR, ILR,  medium,  very broad',
             norm = 2.020000E-13,
             mean = 1236.729,
             fwhm = 255.4998,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (100., 2000.),
                       'skew': (1., 1.)},
             fixed = {'norm': True,
                      'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian5',
             norm = 5.474183E-13,
             mean = 1235.996,
             fwhm = 861.4926,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (100., 2000.),
                       'skew': (1., 1.)},
             fixed = {'norm': True,
                      'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian6',
             norm = 3.948799E-12,
             mean = 1235.138,
             fwhm = 3040.59,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (100., 20000.),
                       'skew': (1., 1.)},
             fixed = {'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian7',
             norm = 2.156964E-12,
             mean = 1231.589,
             fwhm = 8133.099,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (500., 50000.),
                       'skew': (1., 1.)},
             fixed = {'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian8',
             norm = 0.536853E-12,
             mean = 1237.643,
             fwhm = 18183.71,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (900., 3000.),
                       'fwhm': (1200., 1.20E5),
                       'skew': (1., 1.)},
             fixed = {'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'N V - 2 NLR, 2 ILR, 2 medium, 1 very broad',
             norm = 1.217935E-14,
             mean = 1259.753,
             fwhm = 255.4998,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (100.,  2000.),
                       'skew': (1., 1.)},
             tied = {'fwhm': lambda m: 1.0 * m['Ly alpha - NLR, ILR,  medium,  very broad'].fwhm},
             fixed = {'norm': True,
                      'mean': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian10',
             norm = 1.217935E-14,
             mean = 1263.803,
             fwhm = 255.4998,
             bounds = {'norm': (0., 1.00E-11),
                       'mean': (1150., 2000.),
                       'fwhm': (30.,  12000.),
                       'skew': (1., 1.)},
             tied = {'norm': lambda m: 1.0      * m['N V - 2 NLR, 2 ILR, 2 medium, 1 very broad'].norm,
                     'mean': lambda m: 1.003215 * m['N V - 2 NLR, 2 ILR, 2 medium, 1 very broad'].mean,
                     'fwhm': lambda m: 1.0      * m['N V - 2 NLR, 2 ILR, 2 medium, 1 very broad'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian11',
             norm = 6.219548E-15,
             mean = 1259.533,
             fwhm = 861.4926,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (100.,  2000.),
                       'skew': (1., 1.)},
             tied = {'mean': lambda m: 1.019043 * m['gaussian5'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian5'].fwhm},
             fixed = {'norm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian12',
             norm = 6.219548E-15,
             mean = 1263.582,
             fwhm = 861.4926,
             bounds = {'norm': (0., 1.00E-11),
                       'mean': (1150., 2000.),
                       'fwhm': (30.,  12000.),
                       'skew': (1., 1.)},
             tied = {'norm': lambda m: 1.0      * m['gaussian11'].norm,
                     'mean': lambda m: 1.003215 * m['gaussian11'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian11'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian13',
             norm = 0.221692E-13,
             mean = 1258.659,
             fwhm = 3040.59,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (500.,  50000.),
                       'skew': (1., 1.)},
             tied = {'mean': lambda m: 1.019043 * m['gaussian6'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian6'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian14',
             norm = 2.221692E-13,
             mean = 1262.705,
             fwhm = 3040.59,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1150., 2000.),
                       'fwhm': (500.,  50000.),
                       'skew': (1., 1.)},
             tied = {'norm': lambda m: 1.0      * m['gaussian13'].norm,
                     'mean': lambda m: 1.003215 * m['gaussian13'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian13'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian15',
             norm = 3.185217E-13,
             mean = 1255.042,
             fwhm = 8133.099,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (500.,  50000.),
                       'skew': (1., 1.)},
             tied = {'mean': lambda m: 1.019043 * m['gaussian7'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian7'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian16',
             norm = 1.185217E-13,
             mean = 1259.077,
             fwhm = 8133.099,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1150., 2000.),
                       'fwhm': (500.,  50000.),
                       'skew': (1., 1.)},
             tied = {'norm': lambda m: 1.0      * m['gaussian15'].norm,
                     'mean': lambda m: 1.003215 * m['gaussian15'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian15'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'gaussian17',
             norm = 2.287816E-13,
             mean = 1263.24,
             fwhm = 18183.71,
             bounds = {'norm': (0., 1.00E-10),
                       'mean': (1000., 2000.),
                       'fwhm': (1500.,  25000.),
                       'skew': (1., 1.)},
             tied = {'mean': lambda m: 1.020682 * m['gaussian15'].mean,
                     'fwhm': lambda m: 1.0      * m['gaussian15'].fwhm},
             fixed = {'skew': True},
             ) \
+ \
    gaussian(name = 'C III broad absorption',
             norm = -1.000000E-15,
             mean = 1194.,
             fwhm = 3683.102,
             skew = 0.1849483,
             bounds = {'norm': (-1.00E-10, 0.),
                       'mean': (1000., 2000.),
                       'fwhm': (200.,  20000.),
                       'skew': (0.1, 9.)},
             tied = {'fwhm': lambda m: 1.0 * m['Ly alpha broad absorption'].fwhm,
                     'skew': lambda m: 1.0 * m['Ly alpha broad absorption'].skew},
             fixed = {'norm': True,
                      'mean': True},
             ) \
+ \
    gaussian(name = 'Ly alpha broad absorption',
             norm = -5.349239E-13,
             mean = 1235.,
             fwhm = 3683.102,
             skew = 0.1849483,
             bounds = {'norm': (-1.00E-10, 0.),
                       'mean': (1000., 2000.),
                       'fwhm': (200.,  20000.),
                       'skew': (0.1, 9.)},
             fixed = {'mean': True,
                      'fwhm': True,
                      'skew': True},
             ) \
+ \
    gaussian(name = 'N V broad absorption',
             norm = -10.921041E-14,
             mean = 1258.,
             fwhm = 3683.102,
             skew = 0.1849483,
             bounds = {'norm': (-1.00E-10, 0.),
                       'mean': (1000., 2000.),
                       'fwhm': (200.,  20000.),
                       'skew': (0.05, 5.)},
             tied = {'fwhm': lambda m: 1.0 * m['Ly alpha broad absorption'].fwhm,
                     'skew': lambda m: 1.0 * m['Ly alpha broad absorption'].skew},
             fixed = {'mean': True},
             ) \
+ \
    gaussian(name = 'N V broad absorption - 2',
             norm = -8.921041E-14,
             mean = 1262.044,
             fwhm = 3683.102,
             skew = 0.1849483,
             bounds = {'norm': (-1.00E-10, 0.),
                       'mean': (1150., 2000.),
                       'fwhm': (200.,  20000.),
                       'skew': (0.05, 5.)},
             tied = {'norm': lambda m: 1.0      * m['N V broad absorption'].norm,
                     'mean': lambda m: 1.003215 * m['N V broad absorption'].mean,
                     'fwhm': lambda m: 1.0      * m['N V broad absorption'].fwhm,
                     'skew': lambda m: 1.0      * m['N V broad absorption'].skew},
             )
