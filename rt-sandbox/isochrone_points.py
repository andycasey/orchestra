
import numpy as np
from scipy import integrate
from astropy.table import Table
from isochrones.dartmouth import Dartmouth_Isochrone

from smh import radiative_transfer as rt
from smh.linelists import LineList
from smh.photospheres import interpolator as PhotosphereInterpolator


masses = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 5.0]
log10_ages = np.log10([np.arange(0.5, 10, 0.05) * 1e9]).flatten()
metallicities = [-3, -2.5, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]


dartmouth = Dartmouth_Isochrone()

unobservable_parameters = []

stellar_parameters = []
for log10_age in log10_ages:
    for metallicity in metallicities:

        teffs = 10**dartmouth.logTeff(masses, log10_age, metallicity)
        loggs = dartmouth.logg(masses, log10_age, metallicity)
        fehs = metallicity * np.ones(teffs.size)

        stellar_parameters.extend(np.vstack([teffs, loggs, fehs]).T)

        unobservable_parameters.extend(np.vstack([
            masses,
            10**log10_age * np.ones(teffs.size)            
        ]).T)

stellar_parameters = np.array(stellar_parameters)
unobservable_parameters = np.array(unobservable_parameters)

finite = np.all(np.isfinite(stellar_parameters), axis=1)
stellar_parameters = stellar_parameters[finite]
unobservable_parameters = unobservable_parameters[finite]


fig, ax = plt.subplots()
ax.scatter(stellar_parameters[:, 0], stellar_parameters[:, 1],
    c=stellar_parameters[:, 2], edgecolor="none", alpha=0.5)

ax.set_xlim(ax.get_xlim()[::-1])
ax.set_ylim(ax.get_ylim()[::-1])


photosphere_interpolator = PhotosphereInterpolator()
transitions = LineList(
    rows=[(5242.49, 26.000, 3.630, -0.99, 0.0, 0.0, "")], 
    names=("wavelength", "species", "expot", "loggf", "damp_vdw", "dissoc_E", "comments"))


# These are the points to calculate spectra on:
N = len(stellar_parameters)
photospheric_properties = np.nan * np.ones((N, 5, 72), dtype=float)
ew = np.nan * np.ones(N, dtype=float) # [mA]
for i, (teff, logg, feh) in enumerate(stellar_parameters):
    if not (7500 > teff > 4000):
        continue

    photosphere = photosphere_interpolator(teff, logg, feh)
    for j, param in enumerate(photosphere.dtype.names):
        photospheric_properties[i, j, :] = photosphere[param]

    (dispersion, flux, meta), = rt.synthesize(photosphere, transitions)

    ew[i] = integrate.trapz(1.0 - flux, dispersion)

    print(i, N, teff, logg, feh, ew[i])


# Assemble all possible predictors.
t = Table(stellar_parameters, names=("teff", "logg", "feh"))
for i in range(5):
    for j in range(72):
        if i in (2, 3):
            t["p_{}_{}".format(i, j)] = np.log10(photospheric_properties[:, i, j])
        else:
            t["p_{}_{}".format(i, j)] = photospheric_properties[:, i, j]

del t["teff"]

#t["mass"] = unobservable_parameters[:, 0]
#t["age"] = unobservable_parameters[:, 1]
#t["log10_age"] = np.log10(unobservable_parameters[:, 1])


#t["linear_feh"] = 10**t["feh"]
#t["logTeff"] = np.log10(t["teff"])
#t["T^2"] = t["teff"]**2
#t["T^0.5"] = ((t["teff"] - 4000)/2500.0)**0.5
#t["logg^2"] = t["logg"]**2
#t["feh^2"] = t["feh"]**2
#t["inv_loge_teff"] = 1.0/np.log(t["teff"])
#t["g"] = 10**t["logg"]
#t["g_d"] = (t["logg"] >= 3.5).astype(int)


X = np.array([t[k] for k in t.dtype.names]).T
y = ew * 1000


keep = (np.sum(np.isfinite(X), axis=1) > 2) * np.isfinite(y) \
     * (stellar_parameters[:, 0] > 4000) \
     * (stellar_parameters[:, 0] < 6500) \
     * (y < 100)


# only giants 
#keep *= (t["logg"] <= 3.5)

X = X[keep]
y = y[keep]

# Make things isotropic?
for i in range(X.shape[1]):
    offset = np.median(X[:, i])
    scale = np.ptp(X[:, i])
    X[:, i] = (X[:, i] - offset)/scale


from sklearn import linear_model

model = linear_model.LassoCV(cv=20, alphas=10**np.linspace(-3, 3, 10), eps=1e-6)
model.fit(X, y)

fig, ax = plt.subplots()
log_alphas = np.log10(model.alphas_)

ax.plot(log_alphas, model.mse_path_, ":")
ax.plot(log_alphas, model.mse_path_.mean(axis=-1), c='k', lw=2)


fig, ax = plt.subplots()
y2 = model.predict(X)
ax.scatter(y, y2, c=X[:, 1], edgecolor="none", alpha=0.5)

print(np.mean(y2 - y), np.std(y2 - y))

