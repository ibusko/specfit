# Functions for fitting spectral data.

import numpy as np

import astropy.modeling.fitting as fitting
from astropy.modeling import Fittable1DModel, Parameter

from data_objects import SpectrumData

import custom_models
import n5548_models as models


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

    THIS BREAKS WHEN INCLUDING THE EXTINCTION FUNCTION. THE EXTINCTION
    MUST BE MULTIPLIED INTO THE MODEL INSTEAD OF ADDED.

    Parameters
    ----------
    components: list
      list with instances of Fittable1DModel

    Returns
    -------
    An astropy compound model. If no components exist in
    the list, None is returned.

    '''

    # return models.model1

    if len(components) > 0:
        compound_model = components[0]
        for component in components[1:]:
            # composition is for now just additive.
            compound_model += component

        # Set the 'fixed' flag in the compound model parameters.
        # This bug fix is required under astropy 1.0. Under 1.0.1dev
        # this can be removed.
        for component_index in range(len(components)):
            component = components[component_index]
            for parameter_name in component.param_names:
                compound_model_parameter_name = compound_model._param_map_inverse[(component_index,parameter_name)]
                compound_model.fixed[compound_model_parameter_name] = components[component_index].fixed[parameter_name]

        return compound_model
    else:
        return None



tie_call_count = 0

class Tie(object):
    ''' Supplies the tie that links the value of a parameter to
        the value of another parameter.
    '''
    def __init__(self, parameter_name, parent_index, factor):
        self.parameter_name = parameter_name
        self.parent_index = parent_index - 1 # specfit indices are 1-indexed!!!
        self.factor = factor

    # Parameter ties are basically callable objects.
    def __call__(self, *args, **kwargs):
        compound_model = args[0]
        parent_component = compound_model[self.parent_index]
        parent_value = getattr(parent_component, self.parameter_name).value

        tied_value = parent_value * self.factor

        global tie_call_count
        tie_call_count += 1

        return tied_value


def _set_special_attributes(component, fixed, ties):
    for parameter_name in fixed:
        parameter = getattr(component, parameter_name)
        setattr(parameter, 'fixed', fixed[parameter_name])
    for parameter_name in ties:
        parameter = getattr(component, parameter_name)
        tie_parent = ties[parameter_name][0]
        tie_factor = ties[parameter_name][1]
        setattr(parameter, 'tied', Tie(parameter_name, tie_parent, tie_factor))
    return component


# data used by the _build_component function.
constructors = {
    'gaussian': 'custom_models.gaussian(*pars, name=name)',
    'powerlaw': 'custom_models.powerlaw_2(*pars, name=name)',
    'ccmext': 'custom_models.ccmext(*pars, name=name)'
}
component_types = {
    'powerlaw': ['amplitude', 'x_0', 'alpha'],
    'gaussian': ['norm', 'mean', 'fwhm', 'skew'],
    'ccmext':   ['ebmv', 'rv'],
    'dampabs':  ['dummy', 'dummy', 'dummy'],
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
        fixed = {}
        ties = {}
        for count in range(npar):
            line1 = fp.readline()
            tokens1 = line1.split()
            # parameter attributes
            value = float(tokens1[0])
            fixed_flag = int(tokens1[5])
            fixedp = fixed_flag < 0
            pars.append(value)
            parname = component_types[component_type][count]
            fixed[parname] = fixedp

            # The fixed flag also plays the role of a pointer
            # that points to the parent component, when a
            # parameter is linked to the same parameter in the
            # parent.
            # The step size also serves the dual function of
            # supplying the ratio factor when a parameter is
            # linked to another. Both the ratio and the index
            # of the parent component are passed as a tuple
            # to the function that establishes the ties.
            if fixed_flag > 0:
                ties[parname] = ( fixed_flag, float(tokens1[3]) )

        if component_type in constructors:
            constructor = constructors[component_type]
            component = eval(constructor)
            _set_special_attributes(component, fixed, ties)
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
      list with instances of Fittable1DModel

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
            if component_type_index < len(component_types):
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
    for k in range(len(model1.param_names)):
        p_name = model1.param_names[k]

        value1 = model1._parameters[k]
        value2 = model2._parameters[k]
        diff = (value2 - value1) / value1

        if model2.fixed[p_name]:
            fixed = "F"
        else:
            fixed = " "

        tie = model2.tied[p_name]
        if tie and isinstance(tie, Tie):
            parent_index = str(tie.parent_index)
            tie_factpr = str(tie.factor)
        else:
            parent_index = " "
            tie_factpr = " "

        print("%11s   %9.2e      %9.2e       %7.2f          %s      %s   %s" % (p_name, value1, value2, diff, fixed, parent_index, tie_factpr))


def _print_output(x, fit_result, compound_model, fitter, mask, chisq_in, chisq_out, n_free_par, tie_calls, start_time, end_time):
    # we need much better formatting here, but this
    # should suffice as a rudimentary way to compare
    # results with expected values.
    _print_model(compound_model, "\n\n ********** INPUT MODEL ********** \n\n"),
    _print_model(fit_result, "\n\n ********** FITTED MODEL ********** \n\n"),

    # better formatted output
    print("\n\n\n ********** RELATIVE DIFFERENCE ********** \n")
    print("parameter      input model    output model     difference   fixed    tie (index, factor)\n")
    _print_sidebyside(compound_model, fit_result),

    print("\n\n\n ********** REDUCED CHI SQUARE ********** \n\n")
    print("From input model:  %f" % chisq_in)
    print("From output model: %f" % chisq_out)
    print("Total data points: %d" % len(x))
    print("Data points in wavelength ranges: %d" % np.sum(mask))
    print("Number of free parameters: %d" % n_free_par)
    print("\n")
    print("Number of iterations: %d" % fitter.fit_info['nfev'])
    print("Number of tied parameter calls: %d" % tie_calls)
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
    w = mask.copy()
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
    fit_result = fitter(compound_model, x, y, weights=w, acc=1.E-27, maxiter=6000)
    end_time = time.time()

    print(fitter.fit_info['message'])

    cov = fitter.fit_info['param_cov']
    print(cov)
    for i in range(cov.shape[0]):
        print cov[i,i]

    # Parameter errors will be the square root of the diagonal elements.
    #
    # The order of the diagonal elements is the same as in fit_result.param_names

    fit_errors = {}
    i = 0
    if cov is not None:
        for param_name in fit_result.param_names:
            fixed = fit_result.fixed[param_name]
            tied = fit_result.tied[param_name]
            if not fixed and not tied:
                fit_errors[param_name] = math.sqrt(cov[i,i])
                i += 1

    # Must map compound model's internal parameter names to
    # individual component names.

    param_errors = {}
    for param_name in fit_errors.keys():
        index, target_param_name = fit_result._param_map[param_name]
        component_name = fit_result._submodels_names[index]
        param_errors[(component_name, target_param_name)] = fit_errors[param_name]

    # chi-sq
    if 'fixed' in fit_result.parameter_constraints:
        fix = np.asarray(fit_result.fixed.values())
        n_fixed_parameters = np.sum(np.where(fix, 1, 0))
    else:
        n_fixed_parameters = 0

    if 'tied' in fit_result.parameter_constraints:
        tie = np.asarray(fit_result.tied.values())
        n_tied_parameters = np.sum(np.where(tie, 1, 0))
    else:
        n_tied_parameters = 0

    n_free_par = len(fit_result.parameters) - n_fixed_parameters - n_tied_parameters

    if len(e) > 0:
        chisq_in = chisq(x, y, e, mask, compound_model, n_free_par)
        chisq_out = chisq(x, y, e, mask, fit_result, n_free_par)
    else:
        chisq_in = 0.
        chisq_out = 0.

    # report results
    global tie_call_count
    _print_output(x, fit_result, compound_model, fitter, mask, chisq_in, chisq_out, n_free_par, tie_call_count, start_time, end_time)

