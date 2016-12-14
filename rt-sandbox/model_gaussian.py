"""
Predict line strengths by convolving scaled photospheric conditions with a
single Gaussian.
"""


import numpy as np
from astropy.table import Table
from scipy import (integrate, interpolate, optimize as op)
from sklearn import linear_model
from matplotlib import cm

from smh import radiative_transfer as rt
from smh.linelists import LineList
from smh.photospheres import interpolator as PhotosphereInterpolator


# Get the isochrone points.
isochrones = {
    0.00: "dartmouth-solar-isochrone.iso",
    #-0.25: "dartmouth-m0.25-isochrone.iso",
    #-0.50: "dartmouth-m0.5-isochrone.iso",
    #-0.75: "dartmouth-m0.75-isochrone.iso",
    #-1.00: "dartmouth-m1-isochrone.iso",
    -2.00: "dartmouth-m2-isochrone.iso",
}

stellar_parameters = []
for feh, filename in isochrones.items():

    # Calculate EWs for each point.
    iso_logTeff, loggs = np.loadtxt(filename, usecols=(2, 3), unpack=True)
    teffs = 10**iso_logTeff
    fehs = feh * np.ones(teffs.size)

    stellar_parameters.extend(np.vstack([teffs, loggs, fehs]).T)

stellar_parameters = np.array(stellar_parameters)

fgk = (6500 >= stellar_parameters[:, 0]) * (stellar_parameters[:, 0] >= 4000)
stellar_parameters = stellar_parameters[fgk]


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
for i, (teff, logg, feh) in enumerate(stellar_parameters):

    photosphere = photosphere_interpolator(teff, logg, feh)
    for j, param in enumerate(photosphere.dtype.names):
        photospheric_properties[i, j, :] = photosphere[param][:K]

    (dispersion, flux, meta), = rt.synthesize(photosphere, transitions)

    ew[i] = integrate.trapz(1.0 - flux, dispersion) * 1000.0


# Logarithmically scale any of the photospheric properties?
for column_name in ("RHOX", "P", "XNE", "ABROSS"):
    index = photosphere.dtype.names.index(column_name)
    photospheric_properties[:, index, :] = np.log10(photospheric_properties[:, index, :])

# Place a Gaussian in tau, or in T, or in something else?
gaussian_index = photosphere.dtype.names.index("RHOX")
skip_photospheric_names = ["RHOX"] #"XNE", "T", "RHOX", "P", "ABROSS"]

t = Table(stellar_parameters, names=("teff", "logg", "feh"))

for j, name in enumerate(photosphere.dtype.names):
    
    if name in skip_photospheric_names: continue
    for k in range(K):
        column_name = "{}_{}".format(name, k)
        t[column_name] = photospheric_properties[:, j, k]

for k in range(K):
    column_name = "P*ABROSS_{}".format(k)
    t[column_name] = photospheric_properties[:, 2, k] * photospheric_properties[:, 4, k]

t["feh"] = 10**t["feh"]

# Delete columns we don't want to use.
#del t["feh"]

X = np.array([t[k] for k in t.dtype.names]).T
Y = ew.copy()

# Make things ~isotropic?
for i in range(X.shape[1]):
    if i == 1: continue
    offset = np.median(X[:, i])
    scale = np.ptp(X[:, i])
    if scale > 0:
        X[:, i] = (X[:, i] - offset)/scale

assert np.all(np.isfinite(X)) \
   and np.all(np.isfinite(Y)), "X and Y aren't all finite bro"


# At each mu/sigma, try a linear model.
offset = [i for i, n in enumerate(t.dtype.names) if n.endswith("_0")][0]
def predict_line_strength(theta, cv=20, full_output=False):

    # At each mu/sigma we will do a linear model to see if we can fit it all.

    # theta contains:
    mu, sigma, dmu_dfeh = theta
    if sigma > 0.50 and not full_output:
        return np.nan * np.ones(X.shape[0])


    X_convolved = X.copy()

    for i in range(X.shape[0]):

        mu_i = mu + X[i, 2] * dmu_dfeh
        gaussian_x = photospheric_properties[i, gaussian_index]
        phi = np.exp(-(gaussian_x - mu_i)**2 / (2 * sigma**2))

        # Repeat it so it runs over all photospheric properties.
        R = (X_convolved.shape[1] - offset)/gaussian_x.size
        phi_r = np.repeat(phi, R) if R > 1 else phi

        X_convolved[i, offset:] = X_convolved[i, offset:] * phi_r


    model = linear_model.LassoCV(cv=cv)
    model.fit(X_convolved, Y)

    if not full_output:
        return model.predict(X_convolved)

    return (model.predict(X_convolved), model, X_convolved)


def objective_function(theta):
    residual = predict_line_strength(theta) - Y
    print(theta, np.sum(residual**2))
    return residual

def _scalar_objective_function(theta):
    return np.sum(objective_function(theta)**2)

x0 = [1.6, 0.10]
x0 = [3.7, 1.14]
x0 = [0.43, 0.71]
x0 = [0.52, 0.25]
x0 = [1.0, 0.25]


min_v = photospheric_properties[:, 0, :].min()
max_v = photospheric_properties[:, 0, :].max()
thetas = []
result = []
for xi in np.linspace(min_v, max_v, 70):
    for sigma in (0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.0):
        for dmu in (-2, -1, 0, 1, 2):
            thetas.append([xi, sigma, dmu])
            result.append(_scalar_objective_function(thetas[-1]))

#op_params, _ = op.leastsq(objective_function, x0, maxfev=1000)
op_params = op.fmin(_scalar_objective_function, [1.5, 0.10, 0.5], maxfun=100)


Y_P, model, X_convolved = predict_line_strength(op_params, full_output=True)

for column_name in ("P*ABROSS", ):#photosphere.dtype.names:

    try:
        index = t.dtype.names.index("{}_0".format(column_name))
    except ValueError:
        continue

    fig, ax = plt.subplots()
    ax.plot(model.coef_[index:index + K])
    ax.set_title(column_name)


# Modeling the Gaussian on log10(RHOX) and having log10({P, XNE, ABROSS}) yields
# sum_rms of ~542 with mu = 0.43 and sigma = 0.710

# --> T and RHOX can be removed entirely: all zero coefficients

# Removed P, T, RHOX.
# sum_rms of ~614 with mu = 0.46 and sigma = 0.52

# --> XNE is less smooth than before. Remove XNE and just have P, ABROSS

# Removed XNE, left log10(P) and log10(ABROSS)
# sum_rms of ~578 with mu = 0.52 and sigma = 0.49

# --> P is smooth, ABROSS is smooth-"ish". Remove them both and just have:
#     log10(P) * log10(ABROSS)

# Removed log10(P), log10(ABROSS), and included log10(P)*log10(ABROSS)
# sum_rms of ~395 with mu = -0.12 and sigma = 0.74

# --> P*ABROSS coefficients look like a +emission/-absorption, and the gaussian
#     sigma is so wide that just everything is being included. If the P*ABROSS
#     coefficients also look like +emission/-absorption for other isochrones,
#     then maybe we should model a two-component gaussian for P*ABROSS rather
#     than ABROSS/T?

#     --> Or we should make sure that we are calculating impact on the same tau/T.

# --> Trying same model with [Fe/H] = -2.0 isochrone points
# sum_rms of ~438 with mu =-0.68 and sigma = 1.80

# --> Restrict \sigma to < 0.20 and re-run
# sum_rms of 4531 (good for giants, not dwarfs)

# --> Restrict \sigma to < 0.50 and re-run





fig, ax = plt.subplots()
ax.plot()

#op_params, cov, meta, mesg, code = op.leastsq(
#    objective_function, x0, 
#    maxfev=10000, full_output=True)



# Which mu/sigma had the lowest RMS?