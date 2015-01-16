"""
 Some utility functions for generating and playing with the model grid
"""

import numpy as np

from Starfish.grid_tools import PHOENIXGridInterface, HDF5GridStuffer, HDF5Interface, Instrument


def make_hdf5(filename, model_dir='./', Tvals=None, Z_vals=None, alpha_vals=None, logg_vals=None):
    """

    :param filename: The filename for the HDF5 file
    :param model_dir: The base directory containing Phoenix spectra
    :param Tvals: The temperature values to use in the grid (must be list-like)
    :param Z_vals: The metallicity values to use in the grid (must be list-like)
    :param alpha_vals: The alpha/Fe values to use in the grid (must be list-like)
    :param logg_vals: The logg values to use in the grid (must be list-like)
    :return: HDF5Interface object connected to file you just made
    """
    # Convert the parameters to numpy arrays, if they aren't already
    if not isinstance(Tvals, np.ndarray):
        Tvals = np.array(Tvals)
    if not isinstance(Z_vals, np.ndarray):
        Z_vals = np.array(Z_vals)
    if not isinstance(alpha_vals, np.ndarray):
        alpha_vals = np.array(alpha_vals)
    if not isinstance(logg_vals, np.ndarray):
        logg_vals = np.array(logg_vals)

    # Make a raw grid interface
    pg = PHOENIXGridInterface(base=model_dir, Tvals=Tvals, Z_vals=Z_vals, alpha_vals=alpha_vals, logg_vals=logg_vals)

    # Put it in an HDF5 file
    grid = HDF5GridStuffer(pg, filename)
    grid.process_grid()

    # Return an interface object
    return HDF5Interface(filename)


class IGRINS(Instrument):
    '''IGRINS instrument'''

    def __init__(self, name="IGRINS", FWHM=7.5, wl_range=(13000, 26000)):
        super().__init__(name=name, FWHM=FWHM, wl_range=wl_range)
        # sets the FWHM and wl_range