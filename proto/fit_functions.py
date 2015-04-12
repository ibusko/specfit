# Functions for fitting spectral data.

import numpy as np

import astropy.constants as ac
import astropy.modeling.models as models
import astropy.modeling.fitting as fitting
from astropy.modeling import Fittable1DModel

from data_objects import SpectrumData


class gaussian(models.Gaussian1D):
    def __init__(self, *args):

        flux  = args[0][0]
        mean = args[1][0]
        width_kms = args[2][0]

        # sigma = mean * width_kms / 705951.5
        sigma = mean * width_kms / ac.c.to('km/s').value / 2.354820044
        amp = flux / sigma * 0.3989

        super(gaussian, self).__init__(amp, mean, sigma)

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
        alpha  = args[2][0]

        super(powerlaw, self).__init__(amp, x0, alpha)

        fix_amplitude = args[0][1]
        fix_x_0 = args[1][1]
        fix_alpha = args[2][1]

        self.amplitude.fixed = fix_amplitude
        self.x_0.fixed = fix_x_0
        self.alpha.fixed = fix_alpha


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

        if columns[0] == '#':
            continue

        wa.append(float(columns[0]))
        fl.append(float(columns[1]))
        if len(columns) > 2:
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

    else:
        fmask = np.ones(len(wave))

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

def chisq(x, y, e, mask, model, nfree):
    chisq = np.power(((y - model(x)) / e), 2)
    chisq = np.sum(chisq * mask)
    npoints = sum(mask)
    return math.sqrt(chisq / (npoints - nfree - 1))


def _print_model(compound_model, heading):
    print(heading)
    if isinstance(compound_model, Fittable1DModel):
        print(compound_model)
    else:
        for model in compound_model:
            print(model)


def _print_sidebyside(model1, model2):
    for k in range(len(model1._param_names)):
        p_name = model1._param_names[k]
        value1 = model1._parameters[k]
        value2 = model2._parameters[k]
        diff = (value2 - value1) / value1
        print(p_name, value1, value2, diff)


def _print_output(x, fit_result, compound_model, fitter, mask, chisq_in, chisq_out, n_free_par, start_time, end_time):
    # we need much better formatting here, but this
    # should suffice as a rudimentary way to compare
    # results with expected values.

    _print_model(compound_model, "\n\n ********** INPUT MODEL ********** \n\n"),

    _print_model(fit_result, "\n\n ********** FITTED MODEL ********** \n\n"),

    _print_sidebyside(compound_model, fit_result),

    print("\n\n\n ********** REDUCED CHI SQUARE ********** \n\n")
    print("From input model:  %f" % chisq_in)
    print("From output model: %f" % chisq_out)
    print("Total data points: %d" % len(x))
    print("Data points in wavelength ranges: %d" % np.sum(mask))
    print("Number of free parameters: %d" % n_free_par)
    print("Number of iterations: %d" % fitter.fit_info['nfev'])
    elapsed_time = end_time - start_time
    print("\nElapsed time in fitter engine: %f sec" % elapsed_time)


def process_data(*args):
    ''' Reads files with spectral data, wavelength ranges,
    and spectral model. Fits a compound model, fits, and
    prints the fitted results.

    Parameters
    ----------
    datadir: str
       file path to where the data is
    datafile: str
       spectral data file name
    regionsfile: str
       wavelength regions file name
    modelfile: str
       file name of spectral model with first guesses

    '''
    datadir = args[0][0]
    datafile = args[0][1]
    regionsfile = args[0][2]
    modelfile = args[0][3]

    rf_regions = None
    if regionsfile:
        rf_regions = datadir + regionsfile

    spectrum, mask = read_file(datadir + datafile, regions=rf_regions)

    x = spectrum.x.data
    y = spectrum.y.data
    e = spectrum.e.data
    w = mask
    if len(e) > 0:
        w /= e
        max_w = np.max(w)
        w /= max_w

    model = read_model(datadir + modelfile)
    if len(model) > 1:
        compound_model = compoundModel(model)
    else:
        compound_model = model[0]

    fitter = fitting.LevMarLSQFitter()

    start_time = time.time()
    fit_result = fitter(compound_model, x, y, weights=w, acc=1.E-7, maxiter=1000)
    end_time = time.time()

    #    print '@@@@@@     line: 274  - ',fitter.fit_info['param_cov']

    # chi-sq
    fix = np.asarray(fit_result.fixed.values())
    n_free_par = sum(np.where(fix, 0, 1))
    if len(e) > 0:
        chisq_in = chisq(x, y, e, mask, compound_model, n_free_par)
        chisq_out = chisq(x, y, e, mask, fit_result, n_free_par)
    else:
        chisq_in = 0.
        chisq_out = 0.

    _print_output(x, fit_result, compound_model, fitter, mask, chisq_in, chisq_out, n_free_par, start_time, end_time)

