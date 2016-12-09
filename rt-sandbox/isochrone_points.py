
import numpy as np
from scipy import integrate, interpolate
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
"""
for i in range(5):
    # ALL in: 1.75
    # 1.95
    # 2.35

    for j in range(72):
        if i in (2, ):
            t["p_{}_{}".format(i, j)] = np.log(photospheric_properties[:, i, j])
        else:
            t["p_{}_{}".format(i, j)] = photospheric_properties[:, i, j]
"""

# THIS WORKS (as the only photospheric predictor) but dunno why:
for j in range(72):
    t["c_{}".format(j)] = \
        (photospheric_properties[:, 4, j]/photospheric_properties[:, 1, j])



#for j in range(72):
#    t["c_{}".format(j)] = \
#        1.0/(photospheric_properties[:, 1, j])

#for j in range(72):
#    t["d_{}".format(j)] = \
#        (photospheric_properties[:, 2, j]/photospheric_properties[:, 3, j])


# 4 / log(3) --> good; 2.55 rms
# 4 / log(2) --> good; 3.5 rms
# 4 / log(1) --> good: 2.50 rms
# 4 / log(0) --> bad: 5.9 rms
# 4 / 0 --> bad: 5.7
# 4 / 1 --> good: 2.53 rms
# 4 / 2 --> bad: 10 rms
# 4 / 3 --> OK 6.5 rms.

# 3 / 0 --> bad: 7.45 rms
# 3 / 1 --> ok: 3.51 rms
# 3 / 2 --> bad 8.5 rms
# log(3) / 2 --> bad: 22 ms
# log(3) / 1 --> bad:
# log(3) / 0 --> OK, very smooth coefficients,.....!
# log(2) / 0 --> OK, very smooth coeffs
# 2 / log(0) --> bad, but zeroy and peaky
# log(4) / 0: bad

# 0 / 1
"""
for i in range(5):
    for j in range(5):
        if i == j: continue
        for k in range(72):

            v = photospheric_properties[:, i, k]/photospheric_properties[:, j, k]
            t["c_{}_{}_{}".format(i, j, k)] = v
                

            #t["d_{}_{}_{}".format(i, j, k)] = \
            #    photospheric_properties[:, j, k]/photospheric_properties[:, i, k]

"""

"""
for j in range(72):
    t["d_{}".format(j)] = \
        np.log(photospheric_properties[:, 3, j]/photospheric_properties[:, 1, j])



for j in range(72):
    # Kappa not absou
    t["b_{}".format(j)] = np.log(
            photospheric_properties[:, 1, j]**(-5.0/2) \
            / photospheric_properties[:, 3, j] \
            * np.exp(0.75 / (8.6173303e-5 * photospheric_properties[:, 1, j])))
"""

"""

#for j in range(72):
#    t["c_{}".format(j)] = photospheric_properties[:, 0, j]


for j in range(72):
    t["c_{}".format(j)] = \
        (photospheric_properties[:, 3, j] * photospheric_properties[:, 1, j]) \
        / (photospheric_properties[:, 2, j] * photospheric_properties[:, 4, j])

"""


"""
# Put the predictions onto common rho sampling.
min_common_tau = np.nanmax(photospheric_properties[:, 0, 0])
max_common_tau = np.nanmin(photospheric_properties[:, 0, -1])

common_tau = np.linspace(min_common_tau, max_common_tau, 72)

pp = {}
for i in range(photospheric_properties.shape[0]):

    x = photospheric_properties[i, 0, :]

    if not np.isfinite(x).all():
        for j in range(1, 5):
            for k in range(common_tau.size):
                label = "p_{}_{}".format(j, k)

                pp.setdefault(label, [])
                pp[label].append(np.nan)

    else:


        for j in range(1, 5):

            if j in (2, 3):
                y = np.log10(photospheric_properties[i, j, :])
            else:
                y = photospheric_properties[i, j, :]

            tck = interpolate.splrep(x, y)
            new_y = interpolate.splev(common_tau, tck, ext=3)

            for k in range(common_tau.size):

                label = "p_{}_{}".format(j, k)

                pp.setdefault(label, [])
                pp[label].append(new_y[k])


for k, v in pp.items():
    t[k] = v

"""

            


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
y = (ew * 1000)


keep = (np.sum(np.isfinite(X), axis=1) > 2) * np.isfinite(y) \
     * (stellar_parameters[:, 0] > 4000) \
     * (stellar_parameters[:, 0] < 6500)

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

model = linear_model.LassoCV(cv=20, eps=1e-6)
model.fit(X, y)

fig, ax = plt.subplots()
log_alphas = np.log10(model.alphas_)

ax.plot(log_alphas, model.mse_path_, ":")
ax.plot(log_alphas, model.mse_path_.mean(axis=-1), c='k', lw=2)


fig, ax = plt.subplots()
y2 = model.predict(X)
ax.scatter(y, y2, c=X[:, 1], edgecolor="none", alpha=0.5)

print(np.mean(y2 - y), np.std(y2 - y))

raise a

c_score = {}
d_score = {}
e_score = {}

for i in range(0, 5):
    for j in range(i, 5):
        if i == j: continue
        idx = t.dtype.names.index("c_{}_{}_0".format(i, j))
        c_score[(i, j)] = np.sum(np.abs(model.coef_[idx:idx+72]))
        d_score[(i, j)] = np.sum(np.abs(np.diff(model.coef_[idx:idx+72])))

        e_score[(i, j)] = d_score[(i, j)]/c_score[(i, j)]

print("C score")
print(c_score.keys()[np.argmax(c_score.values())])
print("D score")
print(d_score.keys()[np.argmin(d_score.values())])


print("E score")
print(e_score.keys()[np.argmin(e_score.values())])



fig, ax = plt.subplots()
x = np.arange(72)
y = [model.coef_[t.dtype.names.index("c_0_1_{}".format(j))] for j in range(72)]
y2 = [model.coef_[t.dtype.names.index("c_3_4_{}".format(j))] for j in range(72)]

ax.plot(x, y, label="d_min")
ax.plot(x, y2, label="d_max")
plt.legend()





raise a

fig, ax = plt.subplots()
for i in range(5):

    try:
        index = t.dtype.names.index("p_{}_0".format(i))
    except:
        print("Skipping", i)
        continue

    x = np.arange(72)
    y = [model.coef_[t.dtype.names.index("p_{}_{}".format(i, j))] \
        for j in range(x.size)]

    ax.plot(x, y, label=str(i))

ax.set_title("P INDEX")
plt.legend()



fig, ax = plt.subplots()
x = np.arange(72)
y = [model.coef_[t.dtype.names.index("c_{}".format(j))] for j in range(72)]
ax.plot(x, y)
ax.set_title("c")




fig, ax = plt.subplots()
x = np.arange(72)
y = [model.coef_[t.dtype.names.index("b_{}".format(j))] for j in range(72)]
ax.plot(x, y)
ax.set_title("b")

raise a

fig, ax = plt.subplots()
x = np.arange(72)
y = [model.coef_[t.dtype.names.index("d_{}".format(j))] for j in range(72)]
ax.plot(x, y)
ax.set_title("d")





