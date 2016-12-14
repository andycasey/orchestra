"""
Fit line strengths as a function of the photospheric conditions.
"""

import numpy as np
from astropy.table import Table
from scipy import (integrate, interpolate, optimize as op)
from sklearn import linear_model
from matplotlib import cm

from smh import radiative_transfer as rt
from smh.linelists import LineList
from smh.photospheres import interpolator as PhotosphereInterpolator



isochrones = {
    0.25: "dartmouth-p0.25-isochrone.iso",
    -0.25: "dartmouth-m0.25-isochrone.iso",
    -0.50: "dartmouth-m0.5-isochrone.iso",
    -0.75: "dartmouth-m0.75-isochrone.iso",
    -1.25: "dartmouth-m1.25-isochrone.iso",
    -1.75: "dartmouth-m1.75-isochrone.iso",
    -2.25: "dartmouth-m2.25-isochrone.iso",
    -1.00: "dartmouth-m1-isochrone.iso",
    -1.50: "dartmouth-m1.5-isochrone.iso",
    -2.00: "dartmouth-m2-isochrone.iso",

}


isochrones = {
    0.00: "dartmouth-solar-isochrone.iso",
    #-0.25: "dartmouth-m0.25-isochrone.iso",
    #-0.50: "dartmouth-m0.5-isochrone.iso",
    #-0.75: "dartmouth-m0.75-isochrone.iso",
    #-1.00: "dartmouth-m1-isochrone.iso",
    #-2.00: "dartmouth-m2-isochrone.iso",

}


taus = {}
phis = {}
models = {}

for feh in sorted(isochrones.keys()):

    filename = isochrones[feh]

    iso_logTeff, iso_logg = np.loadtxt(filename, usecols=(2, 3), unpack=True)
    iso_teff = 10**iso_logTeff

    teff = iso_teff
    logg = iso_logg
    fehs = feh * np.ones(iso_teff.size)


    fgk = (6500 >= teff) * (teff >= 4000)
    teff = teff[fgk]
    logg = logg[fgk]
    fehs = fehs[fgk]

    stellar_parameters = np.vstack([teff, logg, fehs]).T


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
    skip_photospheric_names = []
    #"T", "P", "XNE", "ABROSS"]# "XNE", "ABROSS"]
    logarithmic_columns = ["P", "XNE"]

    #solar = stellar_parameters[:, 2] == 0
    #stellar_parameters = stellar_parameters[solar]
    #photospheric_properties = photospheric_properties[solar]

    t = Table(stellar_parameters, names=("teff", "logg", "feh"))

    for j, name in enumerate(photosphere.dtype.names):
        if j == 0: continue # because skipping the common tau
        if name in skip_photospheric_names: continue

        for k in range(K):
            column_name = "{}_{}".format(name, k)
            t[column_name] = photospheric_properties[:, j, k]

            if name in logarithmic_columns:
                t[column_name] = np.log10(t[column_name])


    for k in range(K):
        column_name = "KAPPA_TEFF_{}".format(k)
        t[column_name] = t["logg"] * (photospheric_properties[:, 4, k] \
                            /  photospheric_properties[:, 1, k])

                       #* np.log(photospheric_properties[:, 0, k])


    # Do we want to remove any predictors?
    del t["teff"]
    del t["feh"]


    X = np.array([t[k] for k in t.dtype.names]).T
    Y = ew.copy()


    # Make things isotropic?
    for i in range(X.shape[1]):
        if i == 1: continue
        offset = np.median(X[:, i])
        scale = np.ptp(X[:, i])
        if scale > 0:
            X[:, i] = (X[:, i] - offset)/scale


    assert np.all(np.isfinite(X)) \
       and np.all(np.isfinite(Y)), "X and Y aren't all finite bro"


    log10_mean_tau = np.mean(np.log10(photospheric_properties[:, 2, :]), axis=0)


    def predict_EW(theta, x):

        #mu_0, A_0, sigma_0 = theta[:3]
        #mu_1, A_1, sigma_1 = theta[3:]
        b, m_logg, mu_0, A_0, sigma_0, mu_1, A_1, sigma_1 = theta

        if not (A_0 > 0) \
        or not (A_1 > 0):
            return np.nan * np.ones(X.shape[0])


        #assert A_0 > 0 and A_1 > 0

        # x contains (per star):
        # - logg
        # - kappa/teff values as a function of tau
        logg = x[:, 0]
        kappa_teff = x[:, 1:]
        
        EW = np.ones(X.shape[0])
        for j in range(X.shape[0]):

            tau = log10_mean_tau

            # Calculate the \phi coefficients.
            phi = \
                + A_0 * np.exp(-(tau - mu_0)**2/(2 * sigma_0**2)) \
                - A_1 * np.exp(-(tau - mu_1)**2/(2 * sigma_1**2))
            EW[j] = m_logg * logg[j] + np.dot(phi, kappa_teff[i].T) + b

        # For some reason the results are *much* better if we use a mean tau for
        # this metallicity, rather than the tau for each star.
        """
        phi = A_0 * np.exp(-(log10_tau - mu_0)**2/(2 * sigma_0**2)) \
                - A_1 * np.exp(-(log10_tau - mu_1)**2/(2 * sigma_1**2))
        EW = m_logg * logg + np.sum(phi * kappa_teff, axis=1) + b
        """

        return EW


    def objective_function(theta, x):
        residual = predict_EW(theta, x) - Y
        print(theta, np.sum(residual**2))
        return residual

    model = linear_model.LassoCV(cv=20)
    model.fit(X, Y)

    
    # Estimate mu, A, etc from model coeffs.
    pos_coeff = model.coef_[1:] > 0

    idx0 = np.argmax(model.coef_[1:][pos_coeff])
    idx1 = np.argmin(model.coef_[1:][~pos_coeff])

    
    mu_0 = log10_mean_tau[pos_coeff][idx0]
    A_0 = model.coef_[1:][pos_coeff][idx0]
    sigma_0 = 0.10

    mu_1 = log10_mean_tau[~pos_coeff][idx1]
    A_1 = np.abs(model.coef_[1:][~pos_coeff][idx1])
    sigma_1 = 0.10
    
    # Get an intercept and m_logg from the model coefficients
    x0 = np.array([
        model.intercept_, model.coef_[0],
        mu_0, A_0, sigma_0,
        mu_1, A_1, sigma_1])


    #x0 = np.array([100.625, -3.0, 1.07, 72, 0.10, 1.27, 67.5, 0.10])

    op_params, cov, meta, mesg, code = op.leastsq(
        objective_function, x0, args=(X, ), 
        maxfev=10000, full_output=True)

    fig, ax = plt.subplots()
    ax.scatter(Y, predict_EW(op_params, X))
    ax.set_title(feh)




    """
    b, m_logg, mu_0, A_0, sigma_0, mu_1, A_1, sigma_1 = op_params
    phi = np.zeros(log10_mean_tau.size, dtype=float) \
        + A_0 * np.exp(-(log10_mean_tau - mu_0)**2/(2 * sigma_0**2)) \
        - A_1 * np.exp(-(log10_mean_tau - mu_1)**2/(2 * sigma_1**2))
    """

    linear_rms = np.std(model.predict(X) - Y)
    phi_rms = np.std(predict_EW(op_params, X) - Y)

    raise a
    
    fig, ax = plt.subplots()
    ax.plot(log10_mean_tau, model.coef_[1:], c='r')
    ax.plot(log10_mean_tau, phi, c='b')
    ax.set_title(feh)


    #assert phi_rms <= linear_rms


    taus[feh] = log10_mean_tau
    phis[feh] = phi
    models[feh] = (model, x0, op_params)

    break


raise a








fig, ax = plt.subplots()

N = len(isochrones)
colors = cm.get_cmap('Set1', N)

for i, feh in enumerate((-2.0, -1.0, 0.0)):
    x = taus[feh]
    y = phis[feh]

    ax.plot(y, c=colors(float(i)/4.0), lw=2, label=feh)




raise a
op_params = np.array([  
     1.01523130e+02,   1.49756284e+01,   1.12458033e+00,
     5.13822596e+02,   3.54542523e-01,   1.13742989e+00,
     5.00125169e+02,   3.68850428e-01,   4.17346307e+01,
     2.71839346e+01,  -3.09919589e-01,  -2.81604616e-01,
    -5.24281182e+02,  -5.03504979e+02,   1.06771553e-01,
     1.18557609e-01])






b, m_logg, mu_0, A_0, sigma_0, mu_1, A_1, sigma_1 = op_params[:8]
db_dfeh, dm_logg_dfeh, dmu_0_dfeh, dmu_1_dfeh, \
dA0_dfeh, dA1_dfeh, dsigma0_dfeh, dsigma1_dfeh = op_params[8:]

fig, ax = plt.subplots()

for i, feh in enumerate(sorted(phis.keys())):

    intercept = b + db_dfeh * feh
    logg_coeff = m_logg + dm_logg_dfeh * feh

    mu0 = mu_0 + dmu_0_dfeh * feh
    mu1 = mu_1 + dmu_1_dfeh * feh

    A0 = A_0 + dA0_dfeh * feh
    A1 = A_1 + dA1_dfeh * feh

    sigma0 = sigma_0 + dsigma0_dfeh * feh
    sigma1 = sigma_1 + dsigma1_dfeh * feh


    # Calculate the \phi coefficients.
    tau = taus[feh]
    phi = \
        + A0 * np.exp(-(tau - mu0)**2/(2 * sigma0**2)) \
        - A1 * np.exp(-(tau - mu1)**2/(2 * sigma1**2))

    ax.plot(tau, phi, c="rgb"[i])

















# Now go hierarchical.

isochrones = {
    +0.25: "dartmouth-p0.25-isochrone.iso",
    +0.00: "dartmouth-solar-isochrone.iso",
    -0.25: "dartmouth-m0.25-isochrone.iso",
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

for k in range(K):
    column_name = "KAPPA_TEFF_{}".format(k)
    t[column_name] = (photospheric_properties[:, 4, k] \
                        /  photospheric_properties[:, 1, k])

                   #* np.log(photospheric_properties[:, 0, k])


# Do we want to remove any predictors?
del t["teff"]

X = np.array([t[k] for k in t.dtype.names]).T
Y = ew.copy()


# Make things isotropic?
for i in range(X.shape[1]):
    if i == 1: continue # skip metallicity.
    offset = np.median(X[:, i])
    scale = np.ptp(X[:, i])
    if scale > 0:
        X[:, i] = (X[:, i] - offset)/scale


assert np.all(np.isfinite(X)) \
   and np.all(np.isfinite(Y)), "X and Y aren't all finite bro"


# Define the log10_mean_tau for different metallicities.
log10_mean_taus = {}
for feh in np.unique(t["feh"]):
    match = (t["feh"] == feh)

    log10_mean_taus[feh] \
        = np.log10(np.mean(photospheric_properties[match, 0, :], axis=0))




log10_mean_tau = np.log10(np.mean(photospheric_properties[:, 0, :], axis=0))

log10_mean_taus = {}
for feh in np.unique(t["feh"]):

    match = (t["feh"] == feh)
    log10_mean_taus[feh] = \
        np.log10(np.mean(photospheric_properties[match, 0, :], axis=0))


def predict_hierarchical_EW(theta, x):

    b, m_logg, mu_0, A_0, sigma_0, mu_1, A_1, sigma_1 = theta[:8]
    db_dfeh, dm_logg_dfeh, dmu_0_dfeh, dmu_1_dfeh, \
    dA0_dfeh, dA1_dfeh, dsigma0_dfeh, dsigma1_dfeh = theta[8:]

    #assert A_0 > 0 and A_1 > 0

    # x contains (per star):
    # - logg
    # - kappa/teff values as a function of tau
    logg = x[:, 0]
    feh = x[:, 1]
    kappa_teff = x[:, 2:]
    

    intercept = b + db_dfeh * feh
    logg_coeff = m_logg + dm_logg_dfeh * feh

    mu0 = mu_0 + dmu_0_dfeh * feh
    mu1 = mu_1 + dmu_1_dfeh * feh

    A0 = A_0 + dA0_dfeh * feh
    A1 = A_1 + dA1_dfeh * feh

    sigma0 = sigma_0 + dsigma0_dfeh * feh
    sigma1 = sigma_1 + dsigma1_dfeh * feh

    EW = np.zeros(x.shape[0], dtype=float)
    
    for i in range(x.shape[0]):

        # Calculate the \phi coefficients.
        tau = log10_mean_taus[feh[i]]
        phi = \
            + A0[i] * np.exp(-(tau - mu0[i])**2/(2 * sigma0[i]**2)) \
            - A1[i] * np.exp(-(tau - mu1[i])**2/(2 * sigma1[i]**2))

        EW[i] = logg_coeff[i] * logg[i] + np.dot(phi, x[i, 2:].T) + intercept[i]


    return EW




def objective_function(theta, x):
    residual = predict_hierarchical_EW(theta, x) - Y
    print(theta, np.sum(residual**2))
    return residual


x0 = np.array([
    # b_0
    100.625, -3.0, 1.07, 72, 0.10, 1.27, 67.5, 0.10, 
    0, 0, 0, 0, 0, 0, 0, 0
    ])

op_params = np.array([  
     1.01523130e+02,   1.49756284e+01,   1.12458033e+00,
     5.13822596e+02,   3.54542523e-01,   1.13742989e+00,
     5.00125169e+02,   3.68850428e-01,   4.17346307e+01,
     2.71839346e+01,  -3.09919589e-01,  -2.81604616e-01,
    -5.24281182e+02,  -5.03504979e+02,   1.06771553e-01,
     1.18557609e-01])

op_params, cov, meta, mesg, code = op.leastsq(
    objective_function, x0, args=(X, ), 
    maxfev=100000, full_output=True)

fig, ax = plt.subplots()
ax.scatter(Y, predict_hierarchical_EW(op_params, X))


# With single tau: 3674.0 and rms 1.05
# Now with mean tau|feh: 3596.99 and rms 1.04


# Predict for low metallicity stars. How well is our linear approximation in
# [Fe/H] for other things?



isochrones = {
    0.25: "dartmouth-p0.25-isochrone.iso",
    0.00: "dartmouth-solar-isochrone.iso",
    -0.25: "dartmouth-m0.25-isochrone.iso",
    -0.50: "dartmouth-m0.5-isochrone.iso",
    -0.75: "dartmouth-m0.75-isochrone.iso",
    -1.00: "dartmouth-m1-isochrone.iso",
    -1.25: "dartmouth-m1.25-isochrone.iso",
    -1.50: "dartmouth-m1.5-isochrone.iso",
    -1.75: "dartmouth-m1.75-isochrone.iso",
    -2.00: "dartmouth-m2-isochrone.iso",
    -2.25: "dartmouth-m2.25-isochrone.iso",
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

for k in range(K):
    column_name = "KAPPA_TEFF_{}".format(k)
    t[column_name] = (photospheric_properties[:, 4, k] \
                        /  photospheric_properties[:, 1, k])

                   #* np.log(photospheric_properties[:, 0, k])


# Do we want to remove any predictors?
del t["teff"]

X = np.array([t[k] for k in t.dtype.names]).T
Y_mp = ew.copy()


# Make things isotropic?
for i in range(X.shape[1]):
    if i == 1: continue # skip metallicity.
    offset = np.median(X[:, i])
    scale = np.ptp(X[:, i])
    if scale > 0:
        X[:, i] = (X[:, i] - offset)/scale



log10_mean_taus = {}
for feh in np.unique(t["feh"]):

    match = (t["feh"] == feh)
    log10_mean_taus[feh] = \
        np.log10(np.mean(photospheric_properties[match, 0, :], axis=0))


def predict_hierarchical_EW(theta, x):

    b, m_logg, mu_0, A_0, sigma_0, mu_1, A_1, sigma_1 = theta[:8]
    db_dfeh, dm_logg_dfeh, dmu_0_dfeh, dmu_1_dfeh, \
    dA0_dfeh, dA1_dfeh, dsigma0_dfeh, dsigma1_dfeh = theta[8:]

    #assert A_0 > 0 and A_1 > 0

    # x contains (per star):
    # - logg
    # - kappa/teff values as a function of tau
    logg = x[:, 0]
    feh = x[:, 1]
    kappa_teff = x[:, 2:]
    

    intercept = b + db_dfeh * feh
    logg_coeff = m_logg + dm_logg_dfeh * feh

    mu0 = mu_0 + dmu_0_dfeh * feh
    mu1 = mu_1 + dmu_1_dfeh * feh

    A0 = A_0 + dA0_dfeh * feh
    A1 = A_1 + dA1_dfeh * feh

    sigma0 = sigma_0 + dsigma0_dfeh * feh
    sigma1 = sigma_1 + dsigma1_dfeh * feh

    EW = np.zeros(x.shape[0], dtype=float)
    
    for i in range(x.shape[0]):

        # Calculate the \phi coefficients.
        tau = log10_mean_taus[feh[i]]
        phi = \
            + A0[i] * np.exp(-(tau - mu0[i])**2/(2 * sigma0[i]**2)) \
            - A1[i] * np.exp(-(tau - mu1[i])**2/(2 * sigma1[i]**2))

        EW[i] = logg_coeff[i] * logg[i] + np.dot(phi, x[i, 2:].T) + intercept[i]


    return EW










# Now predict!













def predict_hierarchical_EW(theta, x):

    # TODO: WRITE DOWN THE LINEAR ALGEBRA

    t_0, t_1 = theta[:2] # intercept hierarchical params
    a_0, a_1 = theta[2:4] # logg_coeff hierarchical params
    mu0_c0, mu0_c1 = theta[4:6] # coefficients for mu_0
    A0_c0, A0_c1 = theta[6:8] # coefficients for amplitude A0
    sigma0_c0, sigma0_c1 = theta[8:10] # coefficients for sigma_0

    mu1_c0, mu1_c1 = theta[10:12] # coefficients for mu_1
    A1_c0, A1_c1 = theta[12:14] # coefficients for amplitude A1
    sigma1_c0, sigma1_c1 = theta[14:16] # coefficients for sigma_1


    # For each star we need to predict the EW.
    # We need to take the metallicity and get the log10_mean_tau for that [Fe/H]
    # and we need to calculate the gaussian positions, widths, amplitudes from
    # their metallicities

    EWs = np.zeros(X.shape[0], dtype=float)
    for i, xi in enumerate(x):

        logg = xi[0]
        feh = xi[1]
        tau = log10_mean_taus[feh]

        intercept = t_0 + t_1 * feh
        logg_coef = a_0 + a_1 * logg

        mu_0 = mu0_c0 + mu0_c1 * feh
        mu_1 = mu1_c0 + mu1_c1 * feh
        A0 = A0_c0 + A0_c1 * feh
        A1 = A1_c0 + A1_c1 * feh
        sigma0 = sigma0_c0 + sigma0_c1 * feh
        sigma1 = sigma1_c0 + sigma1_c1 * feh

        phi = A0 * np.exp(-(tau - mu_0)**2/(2*sigma0**2)) \
            - A1 * np.exp(-(tau - mu_1)**2/(2*sigma1**2))

        kappa_teff = xi[2:]

        EWs[i] = intercept \
               + logg * logg_coef \
               + np.sum(phi * tau)

    return EWs
        
    
def objective_function(theta, x):
    residual = predict_hierarchical_EW(theta, x) - Y
    print(theta, np.sum(residual**2))
    return residual



#x = (t_0, t_1, a_0, a_1, mu0_c0, mu0_c1, A0_c0, A0_c1, sigma0_c0, sigma0_c1,
#        mu1_c0, mu1_c1, A1_c0, A1_c1, sigma1_c0, sigma1_c1)
x0 = np.array([
    100.417, # t_0
    0.0, #17.3348, # t_1 (100.417 - 96.0833)/0.25
    -3.0, # a_0
    0.0,  # a_1 --> THIS IS A GUESS
    1.063, # mu0_c0
    0.0, #0.24, # mu0_c1,
    76.0, # A0_c0
    0.0, #180.0, # A0_c1
    0.10, # sigma0_c0 --> GUESS
    0.0,  # sigma0_c1 --> GUESS
    1.29, # mu1_c0
    0.0, #0.16, # mu1_c1
    66.0, # A1_c0
    0.0, #60, # A1_c1
    0.10, #sigma1_c0
    0.0,  # sigma1_c1
])

op_params, cov, meta, mesg, code = op.leastsq(
    objective_function, op_params, args=(X, ), 
    maxfev=10000, full_output=True)

fig, ax = plt.subplots()
ax.scatter(Y, predict_hierarchical_EW(op_params, X))


raise a












# 670
