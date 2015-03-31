# test script
#
# attempt at reading spectrum
#
import numpy as np
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

    f.close()

    w = np.array(wave)
    f = np.array(flux)
    e = np.array(error)

    spectrum = SpectrumData()

    spectrum.set_x(w, unit='Angstrom')
    spectrum.set_y(f, unit='erg.s^-1.cm^-2.Angstrom^-1')

    return spectrum


if __name__ == "__main__":
    spectrum = read_file("n5548/n5548_mean_g130mb4.asc")

    # now we need to read the model in here .....

