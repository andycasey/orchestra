"""
Fit line strengths as a function of the photospheric conditions.
"""

import numpy as np
from astropy.table import Table
from scipy import (integrate, interpolate, optimize as op)
from sklearn import linear_model

from smh import radiative_transfer as rt
from smh.linelists import LineList
from smh.photospheres import interpolator as PhotosphereInterpolator


isochrones = {
    0.00: "dartmouth-solar-isochrone.iso",
}
teff = []
logg = []
fehs = []

for feh, filename in isochrones.items():
    iso_logTeff, iso_logg = np.loadtxt(filename, usecols=(2, 3), unpack=True)
    iso_teff = 10**iso_logTeff

    teff.extend(iso_teff)
    logg.extend(iso_logg)
    fehs.extend(feh * np.ones(iso_teff.size))

teff = np.array(teff)
logg = np.array(logg)
fehs = np.array(fehs)


fgk = (6500 >= teff) * (teff >= 4000)
teff = teff[fgk]
logg = logg[fgk]
fehs = fehs[fgk]

stellar_parameters = np.vstack([teff, logg, fehs]).T

fig, ax = plt.subplots()
ax.scatter(stellar_parameters[:, 0], stellar_parameters[:, 1],
    c=stellar_parameters[:, 2], edgecolor="none", alpha=0.5)

ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(ax.get_ylim()[::-1])


# Get ready to start synthesising.
photosphere_interpolator = PhotosphereInterpolator()
transitions = LineList(
    rows=[(5242.49, 26.000, 3.630, -0.99, 0.0, 0.0, "")], 
    names=("wavelength", "species", "expot", "loggf", "damp_vdw", "dissoc_E", "comments"))


# These are the points to calculate spectra on:
N = len(stellar_parameters)
K = 72
photospheric_properties = np.nan * np.ones((N, 5, K), dtype=float)
ew = np.nan * np.ones(N, dtype=float) # [mA]
max_depth = np.nan * np.ones(N, dtype=float)
for i, (teff, logg, feh) in enumerate(stellar_parameters):

    photosphere = photosphere_interpolator(teff, logg, feh)
    for j, param in enumerate(photosphere.dtype.names):
        photospheric_properties[i, j, :] = photosphere[param][:K]

    (dispersion, flux, meta), = rt.synthesize(photosphere, transitions)

    ew[i] = integrate.trapz(1.0 - flux, dispersion) * 1000.0
    max_depth[i] = np.nanmin(flux)

    print(i, N, teff, logg, feh, ew[i], max_depth[i])



# Assemble all possible predictors.
skip_photospheric_names = ["T", "P", "XNE", "ABROSS"]# "XNE", "ABROSS"]

t = Table(stellar_parameters, names=("teff", "logg", "feh"))

for j, name in enumerate(photosphere.dtype.names):
    if j == 0: continue # because skipping the common tau
    if name in skip_photospheric_names: continue

    for k in range(K):
        column_name = "{}_{}".format(name, k)
        t[column_name] = common_photospheric_properties[:, j, k]

        if name in logarithmic_columns:
            t[column_name] = np.log10(t[column_name])


for k in range(K):
    column_name = "KAPPA_TEFF_{}".format(k)
    t[column_name] = (photospheric_properties[:, 4, k] \
                        /  photospheric_properties[:, 1, k])

                   #* np.log(photospheric_properties[:, 0, k])


# Do we want to remove any predictors?
del t["teff"]
del t["feh"]

X = np.array([t[k] for k in t.dtype.names]).T
Y = ew.copy()


# Make things isotropic?
for i in range(X.shape[1]):
    offset = np.median(X[:, i])
    scale = np.ptp(X[:, i])
    if scale > 0:
        X[:, i] = (X[:, i] - offset)/scale


assert np.all(np.isfinite(X)) \
   and np.all(np.isfinite(Y)), "X and Y aren't all finite bro"




# Define intercept and logg coefficient for [Fe/H] = 0.0
b = 100.625
m_logg = -3.0

log10_mean_tau = np.log10(np.mean(photospheric_properties[:, 0, :], axis=0))



def predict_EW(theta, x, log10_tau):

    #mu_0, A_0, sigma_0 = theta[:3]
    #mu_1, A_1, sigma_1 = theta[3:]
    b, m_logg, mu_0, A_0, sigma_0, mu_1, A_1, sigma_1 = theta

    #assert A_0 > 0 and A_1 > 0

    # x contains (per star):
    # - logg
    # - kappa/teff values as a function of tau
    logg = x[:, 0]
    kappa_teff = x[:, 1:]
    
    # Calculate the \phi coefficients.
    phi = np.zeros(log10_mean_tau.size, dtype=float) \
        + A_0 * np.exp(-(log10_mean_tau - mu_0)**2/(2 * sigma_0**2)) \
        - A_1 * np.exp(-(log10_mean_tau - mu_1)**2/(2 * sigma_1**2))
    EW = m_logg * logg + np.dot(phi, kappa_teff.T) + b

    # For some reason the results are *much* better if we use a mean tau for
    # this metallicity, rather than the tau for each star.
    """
    phi = A_0 * np.exp(-(log10_tau - mu_0)**2/(2 * sigma_0**2)) \
            - A_1 * np.exp(-(log10_tau - mu_1)**2/(2 * sigma_1**2))
    EW = m_logg * logg + np.sum(phi * kappa_teff, axis=1) + b
    """

    return EW


def objective_function(theta, x, log10_tau):
    residual = predict_EW(theta, x, log10_tau) - Y
    print(theta, np.sum(residual**2))
    return residual


log10_tau = np.log10(photospheric_properties[:, 0, :])

x0 = np.array([100.625, -3.0, 1.07, 72, 0.10, 1.27, 67.5, 0.10])
op_params, cov, meta, mesg, code = op.leastsq(
    objective_function, x0, args=(X, log10_tau), 
    maxfev=10000, full_output=True)

fig, ax = plt.subplots()
ax.scatter(Y, predict_EW(op_params, X, log10_tau))

# 670
