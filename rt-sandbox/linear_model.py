"""
A huge-ass hierarchical linear model for stellar spectra.
"""

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import collections
import itertools
import multiprocessing as mp
import numpy as np
from sklearn import linear_model
from scipy import optimize as op
from astropy.table import Table

# Load in the synthesised line strengths and photospheric properties.
data = Table.read("line-strengths.fits")

THREADS = 6

# Which photospheric parameter should we model the Gaussians on?
K, DEPTH_PARAMETER = (72, "RHOX")



# Logarithm-ify certain columns.
logarithmify_columns = ("RHOX", "P", "XNE")
for c in logarithmify_columns:
    for k in range(K):
        data["{}_{}".format(c, k)] = np.log10(data["{}_{}".format(c, k)])

# Generate new columns based on existing photospheric columns?


# Delete output columns.
Y = data["_ew"]
del data["_ew"], data["_max_depth"]
assert not any([col.startswith("_") for col in data.dtype.names])

X = np.array([data[column_name] for column_name in data.dtype.names])

#only_consider = (data["feh"] == -2.0)
#only_consider = (data["feh"] >= -1) #* (data["feh"] <= 0.0)
#only_consider = (data["feh"] <= -1)
#X = X[:, only_consider]
#Y = Y[only_consider]

# Checks on PTP values.
ptp = np.ptp(X, axis=1)
ptp_idx = np.argmax(ptp)
print("Largest PTP value is {0:.2e} from {1} (index {2})".format(
    ptp[ptp_idx], data.dtype.names[ptp_idx], ptp_idx))

depth_indices = np.array(
    [data.dtype.names.index("{}_{}".format(DEPTH_PARAMETER, i)) for i in range(K)])

# Make everything roughly isotropic.
P, L = X.shape
for p in range(P):
    if p in depth_indices: continue
    offset = np.mean(X[p])
    scale = np.ptp(X[p])
    if scale > 0:
        X[p] = (X[p] - offset)/scale

assert np.all(np.isfinite(X)) and np.all(np.isfinite(Y)), \
       "Are you sure you know what you're doing?"

# Create some arrays just so it is easier to access the values.
teff = X[data.dtype.names.index("teff")]
logg = X[data.dtype.names.index("logg")]
feh = X[data.dtype.names.index("feh")]

Ne_indices = np.array(
    [data.dtype.names.index("XNE_{}".format(i)) for i in range(K)])
abross_indices = np.array(
    [data.dtype.names.index("ABROSS_{}".format(i)) for i in range(K)])
P_indices = np.array(
    [data.dtype.names.index("P_{}".format(i)) for i in range(K)])
T_indices = np.array(
    [data.dtype.names.index("T_{}".format(i)) for i in range(K)])


_theta_name_order = ("intercept", "teff", "logg", "feh",
    "mu_emission_0", "mu_absorption_0", 
    "sigma_0",
    #"sigma_emission_0", "sigma_absorption_0",
    "dsigma_emission_dlogg", "dsigma_absorption_dlogg",
    "dsigma_dlogg",
    "dmu_emission_dlogg", "dmu_absorption_dlogg",
    "amp_emission_abross_0", "amp_absorption_abross_0",
    "amp_emission_P_0", "amp_absorption_P_0",
    "amp_emission_T_0", "amp_absorption_T_0",
    "amp_emission_Ne_0", "amp_absorption_Ne_0",
    "dmu_emission_dfeh", "dmu_absorption_dfeh",
    "damp_emission_abross_dfeh", "damp_absorption_abross_dfeh",
    "damp_emission_P_dfeh", "damp_absorption_P_dfeh",
    "damp_emission_T_dfeh", "damp_absorption_T_dfeh",
    "damp_emission_Ne_dfeh", "damp_absorption_Ne_dfeh")

# wrap the namedtuple so that we can pass kwargs (e.g. some params that
# we may not want to use anymore in the model)
__theta = collections.namedtuple("Theta", _theta_name_order)
def _theta(**kwargs):
    return __theta(**{k: v for k, v in kwargs.items() if k in __theta._fields})


# For packing and unpacking theta.
_pack_theta = lambda theta: _theta(**dict(zip(_theta_name_order, theta)))
_unpack_theta = lambda theta: [getattr(theta, n) for n in _theta_name_order]

mu_lower_bound, mu_upper_bound = (X[depth_indices].min(), X[depth_indices].max())

# Function to model the data.
def predict_line_strength_by_two_component_gaussian(theta, debug=False):
    """
    Predict the strength of a single absorption line.

    :param theta:
        The model coefficients as a list.

    :returns:
        The predicted line strengths for `L` lines.
    """

    # Pack theta into a named tuple so it is easier to work with.
    T = _pack_theta(theta)

    # Calculate the values at each given metallicity.
    mu_emission = T.mu_emission_0 + T.dmu_emission_dfeh * feh 
    mu_absorption = T.mu_absorption_0 + T.dmu_absorption_dfeh * feh

    sigma_emission = sigma_absorption = T.sigma_0
    # NOTE: requiring single sigma for both, and requiring that it doesn't vary
    #       with [Fe/H]

    if not np.all(mu_upper_bound >= mu_emission) \
    or not np.all(mu_emission >= mu_lower_bound):
        return np.nan * np.ones(L)

    amp_emission_abross = np.clip(T.amp_emission_abross_0 + T.damp_emission_abross_dfeh * feh, -np.inf, np.inf)
    amp_absorption_abross = np.clip(T.amp_absorption_abross_0 + T.damp_absorption_abross_dfeh * feh, -np.inf, np.inf)
    amp_emission_P = np.clip(T.amp_emission_P_0 + T.damp_emission_P_dfeh * feh, -np.inf, np.inf)
    amp_absorption_P = np.clip(T.amp_absorption_P_0 + T.damp_absorption_P_dfeh * feh, -np.inf, np.inf)
    amp_emission_T = np.clip(T.amp_emission_T_0 + T.damp_emission_T_dfeh * feh, -np.inf, np.inf)
    amp_absorption_T = np.clip(T.amp_absorption_T_0 + T.damp_absorption_T_dfeh * feh, -np.inf, np.inf)
    amp_emission_Ne = np.clip(T.amp_emission_Ne_0 + T.damp_emission_Ne_dfeh * feh, -np.inf, np.inf)
    amp_absorption_Ne = np.clip(T.amp_absorption_Ne_0 + T.damp_absorption_Ne_dfeh * feh, -np.inf, np.inf)

    A = np.exp(-(X[depth_indices] - mu_absorption)**2 / 2*sigma_absorption**2)
    E = np.exp(-(X[depth_indices] - mu_emission)**2 / 2*sigma_emission**2)


    EW = np.sum(T.intercept + T.teff * teff + T.logg * logg + T.feh * feh \
       + (A * amp_absorption_abross - E * amp_emission_abross) * X[abross_indices] \
       + (A * amp_absorption_P - E * amp_emission_P) * (X[P_indices] * X[T_indices] * X[Ne_indices]) \
       + (A * amp_absorption_T - E * amp_emission_T) * X[T_indices] \
       + (A * amp_absorption_Ne - E * amp_absorption_Ne) * X[Ne_indices], axis=0)

    if debug:
        raise a
    # TODO: Don't allow negative line strengths?
    return EW


predict_line_strength = predict_line_strength_by_two_component_gaussian

def _objective_function(*args):
    residual = predict_line_strength(*args) - Y
    print(args, np.std(residual))
    return residual

def _scalar_objective_function(*args):
    return np.sum(_objective_function(*args)**2)

x0 = {
    "intercept": np.mean(Y),
    "teff": 0,
    "logg": 0,
    "feh": 0,
    "dsigma_dlogg": 0.0,
    "mu_emission_0": 1.24,
    "mu_absorption_0": 1.40,
    "dsigma_emission_dlogg": 0.0,
    "dsigma_absorption_dlogg": 0.0,
    "dmu_emission_dlogg": 0.0,
    "dmu_absorption_dlogg": 0.0,
    "sigma_0": 0.10,
    "sigma_absorption_0": 0.10,
    "sigma_emission_0": 0.10,
    "amp_emission_abross_0": 1,
    "amp_absorption_abross_0": 1,
    "amp_emission_P_0": 1,
    "amp_absorption_P_0": 1,
    "amp_emission_T_0": 1,
    "amp_absorption_T_0": 1,
    "amp_emission_Ne_0": 1,
    "amp_absorption_Ne_0": 1,
    # Derivatives w.r.t. [Fe/H]
    "dmu_emission_dfeh": 0,
    "dmu_absorption_dfeh": 0,
    "damp_emission_abross_dfeh": 0,
    "damp_absorption_abross_dfeh": 0,
    "damp_emission_P_dfeh": 0,
    "damp_absorption_P_dfeh": 0,
    "damp_emission_T_dfeh": 0,
    "damp_absorption_T_dfeh": 0,
    "damp_emission_Ne_dfeh": 0,
    "damp_absorption_Ne_dfeh": 0,
}

x0 = _unpack_theta(_theta(**x0))
_test = predict_line_strength(x0)

op_params, code = op.leastsq(_objective_function, x0, maxfev=10000)

"""
RMS: 2.4735 over all stars.

array([  0.85923125,  -4.50806829,  -0.09120123,   0.56590817,
        -2.39319847,   2.21548088,   0.68206865,   0.        ,
         0.        ,   0.        ,   0.        ,   0.        ,
         0.87448033,   4.10775913,  12.14066168,   5.48513061,
        -1.77023346,   8.95342997,   1.        ,  -3.93395912,
         0.67539547,   0.68012677,   0.18371188,  -2.24786975,
         2.47883326,   4.67083319,  -3.07569771,   1.47286931,
         0.        ,  -4.25641353])
"""

Y_P = predict_line_strength(op_params)

for column_name in ("logg", "feh"):
    fig, ax = plt.subplots()
    scat = ax.scatter(Y, Y_P, c=data[column_name], edgecolor="none")
    cbar = plt.colorbar(scat)
    cbar.set_label(column_name)


