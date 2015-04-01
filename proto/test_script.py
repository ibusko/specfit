# test script
#
# attempt at reading spectrum and building model
#

import numpy as np
import astropy.modeling.models as models

from data_objects import SpectrumData

wave = []
flux = []
error = []

def read_file(file_name):
    ''' Reads ASCII table with three columns (wavelength,
        flux, error).

    Parameters
    ----------
    file_name: str
       the file path and name

    Returns
    -------
      an instance of SpectrumData

    '''
    f = open(file_name, 'r')
    for line in f:
        columns = line.split()

        wave.append(columns[0])
        flux.append(columns[1])
        error.append(columns[2])

    w = np.array(wave)
    f = np.array(flux)
    e = np.array(error)

    spectrum = SpectrumData()

    spectrum.set_x(w, unit='Angstrom')
    spectrum.set_y(f, unit='erg.s^-1.cm^-2.Angstrom^-1')

    return spectrum


def _compoundModel(components):
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
        sum_of_models = components[0]
        for component in components[1:]:
            sum_of_models += component
        return sum_of_models
    else:
        return None

constructors = {
    'gaussian': 'models.Gaussian1D(*pars)',
    'powerlaw': 'models.PowerLaw1D(*pars)'
}
skipped_parameters = {
    'gaussian': 1,
    }

def _build_component(line, fp, component_type):
    tokens = line.split()
    if len(tokens) > 1:
        name = tokens[0]
        npar = int(tokens[1])
        pars = []
        for count in range(npar):
            line1 = fp.readline()
            tokens1 = line1.split()
            # parameter value
            value = float(tokens1[0])
            pars.append(value)

            # need to throw parameters away. astropy functions
            # are not compatible with specfit models.
        try:
            pars = pars[:-skipped_parameters[component_type]]
        except KeyError:
            pass

        try:
            constructor = constructors[component_type]
            component = eval(constructor)
            return name, component
        except KeyError:
            pass

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
            cname, component = _build_component(line, fp, component_types[component_type_index])
            component_type_index += 1
            if component:
                result.append(component)

    return result


if __name__ == "__main__":
    spectrum = read_file("n5548/n5548_mean_g130mb4.asc")

    #    model = read_model("n5548/sfn5548_lyalpha.modified")
    model = read_model("n5548/sfn5548_lyalpha")

