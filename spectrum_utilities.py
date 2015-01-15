"""
Some functions for making and manipulating IGRINS spectra
"""
from astropy.io import fits
from astropy import units as u
import numpy as np

from Starfish.spectrum import DataSpectrum


def _read_igrins_spectrum(filename):
    hdulist = fits.open(filename)
    orders = []
    for ext in hdulist[1:]:
        x = ext.data.field('wavelength') * u.nm.to(u.angstrom)
        y = ext.data.field('flux')
        e = ext.data.field('error')
        c = ext.data.field('continuum')
        orders.append((x, y, e, c))
    return orders


def _trim_igrins_spectrum(orders):
    """
    The Starfish code needs data in a 2D numpy array, which means each order has to be the same length.
    """
    # First, get the minimum order size. We will use that for all of the orders
    order_size = min([o[0].size for o in orders])

    # Trim all the orders by cutting off both ends of the spectrum
    trimmed_orders = []
    for order in orders:
        size = order[0].size
        if size == order_size:
            trimmed_orders.append(order)
            continue
        left_trim = (order_size - size) / 2
        right_trim = order_size - size - left_trim
        print(order_size, size, left_trim, right_trim)
        x = order[0][left_trim:-right_trim]
        y = order[1][left_trim:-right_trim]
        e = order[2][left_trim:-right_trim]
        c = order[3][left_trim:-right_trim]
        trimmed_orders.append((x, y, c, e))
    return trimmed_orders


def _split_by_type(orders):
    """
      We need to have a 2d numpy array for wavelength, another for flux, and another for error
    """
    wl = np.array([o[0] for o in orders])
    fl = np.array([o[1] for o in orders])
    err = np.array([o[2] for o in orders])
    return wl, fl, err


def read(filename):
    """
     Call this function! Reads in the given spectrum, processes as needed, and returns a DataSpectrum object
    :param filename: the filename of the spectrum to read
    :return: DataSpectrum object
    """
    orders = _read_igrins_spectrum(filename)
    trimmed = _trim_igrins_spectrum(orders)
    wl, fl, err = _split_by_type(trimmed)
    return DataSpectrum(wl, fl, err)