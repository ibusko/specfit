# Model for obs2g.dat in new format

from custom_models import gaussian

model1 = gaussian(name = 'gaussian #1',
                  norm = 1.2E-14,
                  mean = 1001,
                  fwhm = 200,
                  skew = 1.1,
                  bounds = {'norm': (0., 5.E-14),
                            'mean': (1000., 1005.),
                            'fwhm': (50., 500.)},
                  fixed = {'norm': False,
                           'mean': False,
                           'fwhm': False,
                           'skew': True},
                  ) \
+ \
         gaussian(name = 'gaussian #2',
                  norm = 1.E-14,
                  mean = 1002.0,
                  fwhm = 200.,
                  skew = 0.9,
                  fixed = {'norm': False,
                           'mean': False,
                           'fwhm': False,
                           'skew': True},
                  bounds = {'norm': (0., 5.E-14)},
                  tied = {'norm': False,
                          'mean': lambda m: 1.0 * m['gaussian #1'].mean,
                          'fwhm': lambda m: 1.2 * m['gaussian #1'].fwhm,
                          'skew': lambda m: 1.0 * m['gaussian #1'].skew},
                  )




# components	2
# 		gaussian
# 		gaussian
# 		gaussian1  4
# 		         1.100000        0. 0. 0. 0.   0
# 		         4.004000        0. 0. 0. 0.   0
# 		         0.005000        0. 0. 0. 0.   0
# 		         1.000000        0. 0. 0. 0.  -1
# 		gaussian2  4
# 		         0.400000        0. 0. 0. 0.   0
# 		         4.012000        0. 0. 0. 0.   0
# 		         0.007000        0. 0. 0. 0.   0
# 		         1.000000        0. 0. 0. 0.  -1
