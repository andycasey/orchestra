import numpy as np
from numpy import sqrt

def filter(wave, center, fwhm):
    # make a triangular filter
    f = 1. - (1./fwhm) * abs(wave - center)
    f[f < 0.0] = 0.0
    return f

def calc_shk(wave, flux, rv):
    # calculate the (Mount Wilson) S_HK value for a given spectrum
    # inputs:
    # wave : wavelength in Angstrom
    # flux : observed spectral flux at wave
    # rv : stellar RV in km/s; doesn't have to be precise
    # returns:
    # S_MW : S_HK value on the Mount Wilson scale
    # err_S_MW : photon noise error estimate on S_MW
    
    wave = wave/(1. + rv/299792.5)  # shift to star's rest frame
    ind = (wave > 3850.0) & (wave < 4050.0)
    wave = wave[ind]
    flux = flux[ind]
    
    #construct bandpass "filters":
    HKfwhm = 1.09 # width of core filters
    VRfwhm = 20.0 # width of continuum filters
    Kfilter = filter(wave,3933.664,HKfwhm)
    Hfilter = filter(wave,3968.470,HKfwhm)
    Vfilter = filter(wave,3901.070,VRfwhm)
    Rfilter = filter(wave,4001.070,VRfwhm)

    #get flux & photon error in each bandpass:
    Kflux = sum(flux * Kfilter)
    Kerr = sqrt(sum((sqrt(flux) * Kfilter)**2))
    Hflux = sum(flux * Hfilter)
    Herr = sqrt(sum((sqrt(flux) * Hfilter)**2))
    Vflux = sum(flux * Vfilter)
    Verr = sqrt(sum((sqrt(flux) * Vfilter)**2))
    Rflux = sum(flux * Rfilter)
    Rerr = sqrt(sum((sqrt(flux) * Rfilter)**2))

    S_HARPS = (Hflux/HKfwhm + Kflux/HKfwhm)/(Rflux/VRfwhm + Vflux/VRfwhm)
    err_S_HARPS = S_HARPS * sqrt(((Kerr/HKfwhm)**2 + (Herr/HKfwhm)**2)/(Hflux/HKfwhm + Kflux/HKfwhm)**2 + 
                ((Rerr/VRfwhm)**2 + (Verr/VRfwhm)**2)/(Rflux/VRfwhm + Vflux/VRfwhm)**2)
    S_MW = 1.111*S_HARPS + 0.0153  #conversion to Mount Wilson value from Lovis+2011
    err_S_MW = 1.111*err_S_HARPS
    
    return S_MW, err_S_MW
    