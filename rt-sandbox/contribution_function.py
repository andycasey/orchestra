

import numpy as np
from astropy.table import Table
from scipy import integrate, interpolate
from sklearn import linear_model

from smh import radiative_transfer as rt
from smh.linelists import LineList
from smh.photospheres import interpolator as PhotosphereInterpolator


isochrones = {
    0.25: "dartmouth-p0.25-isochrone.iso",
    0.00: "dartmouth-solar-isochrone.iso",
    -0.25: "dartmouth-m0.25-isochrone.iso",
    #-0.50: "dartmouth-m0.5-isochrone.iso",
    #-0.75: "dartmouth-m0.75-isochrone.iso",
    #-1.00: "dartmouth-m1-isochrone.iso",
    #-1.25: "dartmouth-m1.25-isochrone.iso",
    #-1.50: "dartmouth-m1.5-isochrone.iso",
    #-1.75: "dartmouth-m1.75-isochrone.iso",
    #-2.00: "dartmouth-m2-isochrone.iso",
    #-2.25: "dartmouth-m2.25-isochrone.iso",
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
K = 62
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


# Define a common tau scale which we will interpolate on to.
logarithmic_columns = []
min_common_tau = np.nanmax(photospheric_properties[:, 0, 0])
max_common_tau = np.nanmin(photospheric_properties[:, 0, -1])
common_tau = np.linspace(min_common_tau, max_common_tau, K)

common_photospheric_properties = np.nan * np.ones_like(photospheric_properties)
for i in range(N):

    x = photospheric_properties[i, 0, :]
    assert np.isfinite(x).all(), "err, what?"

    # Just for my autisms
    common_photospheric_properties[i, 0, :] = common_tau

    for j in range(1, 5):

        column = photosphere.dtype.names[j]
        y = photospheric_properties[i, j, :]

        if column in logarithmic_columns:
            y = np.log10(y)

        tck = interpolate.splrep(x, y)
        new_y = interpolate.splev(common_tau, tck, ext=3)

        if column in logarithmic_columns:
            new_y = 10**new_y

        common_photospheric_properties[i, j, :] = new_y


logg_coeffs = {}
intercepts = {}

common_fig, common_ax = plt.subplots()

for feh in np.sort(np.unique(stellar_parameters[:, 2])):

    # Assemble all possible predictors.
    skip_photospheric_names = ["T", "P", "XNE", "ABROSS"]# "XNE", "ABROSS"]

    t = Table(stellar_parameters, names=("teff", "logg", "feh"))

    for j, name in enumerate(photosphere.dtype.names):
        if j == 0: continue # because skipping the common tau
        if name in skip_photospheric_names: continue

        for k in range(common_tau.size):
            column_name = "{}_{}".format(name, k)
            t[column_name] = common_photospheric_properties[:, j, k]

            if name in logarithmic_columns:
                t[column_name] = np.log10(t[column_name])


    for k in range(common_tau.size):
        column_name = "4/1_{}".format(k)
        t[column_name] = (photospheric_properties[:, 4, k] \
                            /  photospheric_properties[:, 1, k]) \
                       * np.log(photospheric_properties[:, 0, k])



    """
    for k in range(common_tau.size):
        column_name = "{}_{}".format("PT/XNE", k)
        t[column_name] = (common_photospheric_properties[:, 1, k] \
                       /  common_photospheric_properties[:, 2, k]) \
                       *  np.log(common_photospheric_properties[:, 3, k])
    """

    # Amplify signals
    #for k in range(common_tau.size):
    #    column_name = "{}_{}".format("SIGNAL", k)
    #    t[column_name] = \
    #        common_photospheric_properties[:, 1, k] \
    #      * common_photospheric_properties[:, 3, k] \
    #      * common_photospheric_properties[:, 4, k] \
    #      / common_photospheric_properties[:, 2, k]


    #for k in range(common_tau.size):
    #    column_name = "1/LOG(T)_{}".format(k)
    #    t[column_name] = 1.0/np.log(common_photospheric_properties[:, 1, k])

    #for k in range(common_tau.size):
    #    column_name = "LOG(P)_{}".format(k)
    #    t[column_name] = np.log(common_photospheric_properties[:, 2, k])


    #for k in range(common_tau.size):
    #    column_name = "LOG(P)*LOG(ABROSS)*LOG(T)_{}".format(k)
    #    t[column_name] = (np.log(common_photospheric_properties[:, 2, k]) \
    #                   *  np.log(common_photospheric_properties[:, 1, k]) \
    #                   *  np.log(common_photospheric_properties[:, 3, k]))\
    #                   * np.log(common_photospheric_properties[:, 4, k])


    """
    # Add in a psuedo-Kappa
    name = "KAPPA"
    for k in range(common_tau.size):
        column_name = "{}_{}".format(name, k)
        t[column_name] = np.log(
                common_photospheric_properties[:, 1, j]**(-5.0/2) \
                / common_photospheric_properties[:, 3, j] \
                * np.exp(0.75 / (8.6173303e-5 * common_photospheric_properties[:, 1, j])))

    """

    # Do we want to remove any predictors?
    del t["teff"]

    X = np.array([t[k] for k in t.dtype.names]).T
    Y = ew.copy()


    assert np.all(np.isfinite(X)) \
       and np.all(np.isfinite(y)), "X and Y aren't all finite bro"

    keep = (t["feh"] == feh)
    X = X[keep, :]
    Y = Y[keep]


    # Make things isotropic?
    for i in range(X.shape[1]):
        offset = np.median(X[:, i])
        scale = np.ptp(X[:, i])
        if scale > 0:
            X[:, i] = (X[:, i] - offset)/scale


    model = linear_model.LassoCV(cv=20)
    model.fit(X, Y)

    fig, ax = plt.subplots()
    log_alphas = np.log10(model.alphas_)

    ax.plot(log_alphas, model.mse_path_, ":")
    ax.plot(log_alphas, model.mse_path_.mean(axis=-1), c='k', lw=2)
    ax.set_title(feh)

    fig, ax = plt.subplots()
    Y_P = model.predict(X)
    ax.scatter(Y, Y_P, c=X[:, 1], edgecolor="none", alpha=0.5)
    ax.set_title(feh)
    print(np.mean(Y_P - Y), np.std(Y_P - Y))


    # Look at the coefficients.
    fig, ax = plt.subplots()
    for name in photosphere.dtype.names:
        if "{}_0".format(name) in t.dtype.names:

            y = [model.coef_[t.dtype.names.index("{}_{}".format(name, j))] \
                    for j in range(common_tau.size)]

            ax.plot(common_tau, y, label=name)

        else:
            print("Skipping {0} because {0}_0 not in table".format(name))

    # Show coefficients of custom relations we tried
    for name in t.dtype.names:
        if name.endswith("_0") and name[:-2] not in photosphere.dtype.names:
            y = [model.coef_[t.dtype.names.index("{}_{}".format(name[:-2], j))] \
                    for j in range(common_tau.size)]
            ax.plot(common_tau, y, label=name[:-2])


    plt.legend()
    ax.set_title(feh)


    fig, ax = plt.subplots()
    name = "4/1"

    y = np.array(
        [model.coef_[t.dtype.names.index("{}_{}".format(name, j))] \
             * t["{}_{}".format(name, j)] for j in range(common_tau.size)])

    for i in range(y.shape[1]):
        ax.plot(common_tau, y[:, i], c="k", alpha=0.01)

    ax.set_title(feh)


    c = np.array(
        [model.coef_[t.dtype.names.index("{}_{}".format(name, j))] \
            for j in range(common_tau.size)])

    common_ax.plot(
        np.log10(np.mean(photospheric_properties[:, 0, :], axis=0)), c, label=feh)


    logg_coeffs[feh] = model.coef_[0]
    intercepts[feh] = model.intercept_



fig, ax = plt.subplots()
ax.scatter(logg_coeffs.keys(), logg_coeffs.values())


fig, ax = plt.subplots()
ax.scatter(intercepts.keys(), intercepts.values())














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





