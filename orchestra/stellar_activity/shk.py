
"""
Calculate the S_HK stellar activity index on the Mount Wilson scale.
"""

__authors__ = [
    "Megan Bedell <mbedell@oddjob.uchicago.edu>",
    # With minor docstring/speed modifications from:
    "Andy Casey <arc@ast.cam.ac.uk>"
]


import numpy as np

def _triangular_filter(wave, center, fwhm):
    """
    Construct a triangular filter about some central wavelength.

    :param wave:
        Wavelength array in Angstroms.

    :param center:
        Central wavelength for the filter.

    :param fwhm:
        The full-width half-maximum of the filter.

    :returns:
        An array the same length of `wave` with entries between 0.0 and 1.0.
    """

    return np.clip(1.0 - (1./fwhm) * np.abs(wave - center), 0.0, 1.0)


def shk_index(wave, flux, rv):
    """
    Calculate the (Mount Wilson) S_HK value for a given spectrum.

    :param wave:
        Wavelength array in Angstroms.

    :param flux:
        Observed flux at a given wavelength.

    :param rv:
        Stellar radial velocity in km/s. Value doesn't have to be precise.

    :returns:
        A two-length tuple containing the S_HK value on the Mount Wilson scale,
        and the photon noise error estimate on that S_HK value.
    """

    wave = wave/(1. + rv/299792.458)  # shift to star's rest frame
    ind = (wave > 3850.0) & (wave < 4050.0)
    # Exclude negative or non-finite fluxes in our windowed region
    ind[ind] *= (flux[ind] > 0) & np.isfinite(flux[ind])
    wave = wave[ind]
    flux = flux[ind]
    
    # Construct bandpass "filters":
    HKfwhm = 1.09 # Width of core filters
    VRfwhm = 20.0 # Width of continuum filters
    Kfilter = _triangular_filter(wave, 3933.664, HKfwhm)
    Hfilter = _triangular_filter(wave, 3968.470, HKfwhm)
    Vfilter = _triangular_filter(wave, 3901.070, VRfwhm)
    Rfilter = _triangular_filter(wave, 4001.070, VRfwhm)

    # Get flux & photon error in each bandpass:
    sqrt_flux = flux**0.5
    Kflux = np.sum(flux * Kfilter)
    Kerr = (np.sum((sqrt_flux * Kfilter)**2))**0.5
    Hflux = np.sum(flux * Hfilter)
    Herr = (np.sum((sqrt_flux * Hfilter)**2))**0.5
    Vflux = np.sum(flux * Vfilter)
    Verr = (np.sum((sqrt_flux * Vfilter)**2))**0.5
    Rflux = np.sum(flux * Rfilter)
    Rerr = (np.sum((sqrt_flux * Rfilter)**2))**0.5

    S_HARPS = (Hflux/HKfwhm + Kflux/HKfwhm)/(Rflux/VRfwhm + Vflux/VRfwhm)
    err_S_HARPS = S_HARPS * (
        ((Kerr/HKfwhm)**2 + (Herr/HKfwhm)**2)/(Hflux/HKfwhm + Kflux/HKfwhm)**2 \
      + ((Rerr/VRfwhm)**2 + (Verr/VRfwhm)**2)/(Rflux/VRfwhm + Vflux/VRfwhm)**2)**0.5

    # Conversion to Mount Wilson value from Lovis et al. (2011)
    S_MW = 1.111*S_HARPS + 0.0153
    err_S_MW = 1.111*err_S_HARPS
    
    return (S_MW, err_S_MW)
    