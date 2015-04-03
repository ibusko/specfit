# test script
#
# attempt at reading spectrum and building model
#

import numpy as np

import astropy.modeling.models as models
import astropy.modeling.fitting as fitting

from data_objects import SpectrumData


class gaussian(models.Gaussian1D):
    def __init__(self, *args):

        amp  = args[0][0]
        mean = args[1][0]
        stdd = args[2][0]

        super(gaussian, self).__init__(amp, mean, stdd)

        fix_amplitude = args[0][1]
        fix_mean = args[1][1]
        fix_stdd = args[2][1]

        self.amplitude.fixed = fix_amplitude
        self.mean.fixed = fix_mean
        self.stddev.fixed = fix_stdd


class powerlaw(models.PowerLaw1D):
    def __init__(self, *args):

        amp = args[0][0]
        x0  = args[1][0]

        super(powerlaw, self).__init__(amp, x0)

        fix_amplitude = args[0][1]
        fix_x_0 = args[1][1]

        self.amplitude.fixed = fix_amplitude
        self.x_0.fixed = fix_x_0


def read_file(file_name, regions=None):
    ''' Reads ASCII table with three columns (wavelength,
        flux, error). Can also read a table with wavelength
        regions that define valid data.

    Parameters
    ----------
    file_name: str
       the file path and name
    regions: str, optional
       an ASCII file with wavelength ranges that define
       valid data, one range per row

    Returns
    -------
      an instance of SpectrumData

    '''
    wa = []
    fl = []
    er = []

    f = open(file_name, 'r')
    for line in f:
        columns = line.split()

        wa.append(float(columns[0]))
        fl.append(float(columns[1]))
        er.append(float(columns[2]))

    wave = np.array(wa)
    flux = np.array(fl)
    error = np.array(er)

    spectrum = SpectrumData()

    spectrum.set_x(wave, unit='Angstrom')
    spectrum.set_y(flux, unit='erg.s^-1.cm^-2.Angstrom^-1')
    spectrum.set_e(error, unit='erg.s^-1.cm^-2.Angstrom^-1')

    # Note that SpectrumData does not use masked arrays.
    if regions:
        mask = np.zeros(len(wave))

        f = open(regions, 'r')
        for line in f:
            region = line.split()
            index1 = np.where(wave > float(region[0]))
            index2 = np.where(wave <  float(region[1]))

            mask1 = np.zeros(len(wave))
            mask1[index1] = 1
            mask2 = np.zeros(len(wave))
            mask2[index2] = 1

            mask3 = np.logical_and(mask1, mask2)
            mask = np.logical_or(mask, mask3)

        # Make the boolean mask into a float array so it can
        # be used directly as a weight array by the fitter.
        fmask = np.where(mask, 1.0, 0.0)

        # Can't set a mask in a NDData object. The mask setting causes
        # the shape of the ndarray inside the NDData object to change to
        # an (apparently) arbitrary value. For now, we pass the mask
        # as an independent entity so we can continue to progress on
        # the main task.
        #
        # spectrum.mask = fmask

    return spectrum, fmask


def compoundModel(components):
    ''' Builds the compound model for the active list of components.

    For now, a simple additive composition is used. This will have
    to be modified eventually so as to enable any kind of allowable
    composition.

    Parameters
    ----------
    components: list
      list with instances of Fittable1DModel

    Returns
    -------
    An astropy compound model. If no components exist in
    the list, None is returned.

    '''
    if len(components) > 0:
        compound_model = components[0]
        for component in components[1:]:
            # composition is for now just additive.
            compound_model += component

        # set the 'fixed' flag in the compound model parameters.
        for component_index in range(len(components)):
            component = components[component_index]
            for parameter_name in component.param_names:
                compound_model_parameter_name = compound_model._param_map_inverse[(component_index,parameter_name)]
                compound_model.fixed[compound_model_parameter_name] = components[component_index].fixed[parameter_name]

        return compound_model
    else:
        return None


# data used by the _build_component function.
constructors = {
    'gaussian': 'gaussian(*pars)',
    'powerlaw': 'powerlaw(*pars)'
}
discarded_parameters = {
    'gaussian': 1,
}
first = True

def _build_component(line, fp, component_type):
    global first
    if first:
        first = False
        line1 = fp.readline()
        tokens = line1.split()
    else:
        tokens = line.split()

    if len(tokens) > 1:
        name = tokens[0]
        npar = int(tokens[1])
        pars = []
        for count in range(npar):
            line1 = fp.readline()
            tokens1 = line1.split()
            # parameter attributes
            value = float(tokens1[0])
            fixed = int(tokens1[5]) < 0
            pars.append( (value, fixed) )

        # need to throw parameters away. astropy functions
        # are not directly compatible with specfit models.
        if component_type in discarded_parameters:
            pars = pars[:-discarded_parameters[component_type]]

        if component_type in constructors:
            constructor = constructors[component_type]
            component = eval(constructor)
            return name, component

    return None, None


def read_model(file_name):
    ''' Reads ASCII file with spectral model, in the
        format accepted by the specfit IRAF task.

    Parameters
    ----------
    file_name: str
       the file path and name

    Returns
    -------
      list with instances of Fittable1DModel and Polynomial1DModel

    '''
    n_components = 0
    component_type_index = 0
    component_types = []
    stop_char = '_'
    result = []
    fp = open(file_name, 'r')

    while True:
        line = fp.readline()
        if not line or line.startswith(stop_char):
            break
        tokens = line.split()

        # read component types
        if len(tokens) > 0 and tokens[0] == 'components':
            stop_char = '#'
            n_components = int(tokens[1])
            for count in range(n_components):
                line = fp.readline()
                component_types.append(line.split()[0])

        # read each component
        if n_components > 0:
            cname, component = _build_component(line, fp, component_types[component_type_index])
            component_type_index += 1
            if component:
                result.append(component)

    return result


import math
import time

def _chisq(x, y, e, mask, model, nfree):
    chisq = np.power(((y - model(x)) / e), 2)
    chisq = np.sum(chisq * mask)
    npoints = sum(mask)
    return math.sqrt(chisq / (npoints - nfree))

if __name__ == "__main__":

    datadir = "../data/n5548/"

    spectrum, mask = read_file(datadir + "n5548_mean_g130mb4.asc", regions=datadir + "n5548_lyalpha_sample.dat")
    x = spectrum.x.data
    y = spectrum.y.data
    e = spectrum.e.data

    w = mask / e

    model = read_model(datadir + "sfn5548_lyalpha2")
    compound_model = compoundModel(model)

    fitter = fitting.LevMarLSQFitter()

    start_time = time.time()
    fit_result = fitter(compound_model, x, y, weights=w, acc=1.E-6, maxiter=1000)
    end_time = time.time()

#    print '@@@@@@     line: 274  - ',fitter.fit_info['param_cov']

    # chi-squared
    fix = np.asarray(fit_result.fixed.values())
    free = np.where(fix, 0, 1)
    n_free = sum(free)

    chisq_in  = _chisq(x, y, e, mask, compound_model, n_free)
    chisq_out = _chisq(x, y, e, mask, fit_result, n_free)

    # we need much better formatting here, but this
    # should suffice as a rudimentary way to compare
    # results with expected values.
    print("\n\n ********** INPUT MODEL ********** \n\n")
    for model in compound_model:
        print(model)
    print("\n\n\n ********** FITTED MODEL ********** \n\n")
    for model in fit_result:
        print(model)
    print("\n\n\n ********** REDUCED CHI SQUARE ********** \n\n")
    print("From input model:  %f" % chisq_in)
    print("From output model: %f" % chisq_out)
    print("Total data points: %d" % len(x))
    print("Data points in wavelength ranges: %d" % np.sum(mask))
    print("Number of free parameters: %d" % n_free)

    elapsed_time = end_time - start_time
    print("\nElapsed time in fitter engine: %f sec" % elapsed_time)




