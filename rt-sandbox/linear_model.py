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

# Which photospheric parameter should we model the Gaussians on?
K, DEPTH_PARAMETER = (72, "RHOX")

# Logarithm-ify certain columns.
logarithmify_columns = ("RHOX", "XNE", "P", )
for c in logarithmify_columns:
    for k in range(K):
        data["{}_{}".format(c, k)] = np.log10(data["{}_{}".format(c, k)])

# Delete output columns.
Y = data["_ew"]
del data["_ew"], data["_max_depth"]
assert not any([col.startswith("_") for col in data.dtype.names])

X = np.array([data[column_name] for column_name in data.dtype.names])

only_consider = np.ones(Y.size, dtype=bool)
X = X[:, only_consider]
Y = Y[only_consider]

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
    "sigma_emission_0", "sigma_absorption_0",
    "dsigma_absorption_dfehdlogg",
    "amp_emission_abross_0", "amp_absorption_abross_0",
    "amp_emission_P_0", "amp_absorption_P_0",
    "amp_emission_T_0", "amp_absorption_T_0",
    "amp_emission_Ne_0", "amp_absorption_Ne_0",
    "dmu_emission_dfeh", "dmu_absorption_dfeh",
    "damp_abross_dfeh", 
    "damp_P_dfeh", 
    "damp_T_dfeh", 
    "damp_Ne_dfeh")

# wrap the namedtuple so that we can pass kwargs (e.g. some params that
# we may not want to use anymore in the model)
_theta = collections.namedtuple("Theta", _theta_name_order)
def Theta(**kwargs):
    # Set the default value to zero.
    kwds = dict(zip(_theta._fields, [0] * len(_theta._fields)))
    kwds.update({k: v for k, v in kwargs.items() if k in kwds})
    return _theta(**kwds)


# For packing and unpacking theta.
_pack_theta = lambda theta: Theta(**dict(zip(_theta_name_order, theta)))
_unpack_theta = lambda theta: [getattr(theta, n) for n in _theta_name_order]

mu_lower_bound, mu_upper_bound = (X[depth_indices].min(), X[depth_indices].max())

# Function to model the data.
def predict_line_strength_by_two_component_gaussian(theta, full_output=False):
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

    if not np.all(mu_upper_bound >= mu_emission) \
    or not np.all(mu_emission >= mu_lower_bound):
        return np.nan * np.ones(L)

    sigma_emission = T.sigma_emission_0
    sigma_absorption = T.sigma_absorption_0 + T.dsigma_absorption_dfehdlogg * 10**feh * logg
    
    amp_emission_abross = T.amp_emission_abross_0 + T.damp_abross_dfeh * feh
    amp_absorption_abross = T.amp_absorption_abross_0 + T.damp_abross_dfeh * feh
    amp_emission_P = T.amp_emission_P_0 + T.damp_P_dfeh * feh
    amp_absorption_P = T.amp_absorption_P_0 + T.damp_P_dfeh * feh
    amp_emission_T = T.amp_emission_T_0 + T.damp_T_dfeh * feh
    amp_absorption_T = T.amp_absorption_T_0 + T.damp_T_dfeh * feh
    amp_emission_Ne = T.amp_emission_Ne_0 + T.damp_Ne_dfeh * feh
    amp_absorption_Ne = T.amp_absorption_Ne_0 + T.damp_Ne_dfeh * feh

    A = np.exp(-(X[depth_indices] - mu_absorption)**2 / 2*sigma_absorption**2)
    E = np.exp(-(X[depth_indices] - mu_emission)**2 / 2*sigma_emission**2)

    EW = np.sum(T.intercept + T.teff * teff + T.logg * logg + T.feh * feh \
       + (A * amp_absorption_abross - E * amp_emission_abross) * X[abross_indices] \
       + (A * amp_absorption_P      - E * amp_emission_P) * X[P_indices] \
       + (A * amp_absorption_T      - E * amp_emission_T) * X[T_indices] \
       + (A * amp_absorption_Ne     - E * amp_emission_Ne) * X[Ne_indices], axis=0)

    # TODO: Don't allow negative EWs?

    return (EW, X[depth_indices], A, E) if full_output else EW



predict_line_strength = predict_line_strength_by_two_component_gaussian

def _objective_function(*args):
    residual = predict_line_strength(*args) - Y
    print(args, np.std(residual))
    return residual

def _scalar_objective_function(*args):
    return np.sum(_objective_function(*args)**2)


# Generate initial value.
x0 = _unpack_theta(Theta(
    intercept=0.85972998553986724,
    teff=0.017460341505266527,
    logg=-4.7019731537599894,
    feh=1.0584179693653297,
    mu_emission_0=1.0709724696264904,
    mu_absorption_0=1.1900378387312363,
    sigma_emission_0=0.36229703874314678,
    sigma_absorption_0=0.34797406221173532,
    dsigma_absorption_dfehdlogg=0.04201639903294907,
    amp_emission_abross_0=42.496852825436584,
    amp_absorption_abross_0=44.210684821888307,
    amp_emission_P_0=-15.405794047318949,
    amp_absorption_P_0=-11.752735375086061,
    amp_emission_T_0=72.575016553477411,
    amp_absorption_T_0=70.219395057343903,
    amp_emission_Ne_0=-52.576660952260198,
    amp_absorption_Ne_0=-48.876583107641153,
    dmu_emission_dfeh=-1.0838648723081996,
    dmu_absorption_dfeh=-1.0047589456757466,
    damp_abross_dfeh=-59.361840396035966,
    damp_P_dfeh=-61.768499259040382,
    damp_T_dfeh=-66.199055798272511, 
    damp_Ne_dfeh=86.585760456771794))

#RMS: 1.40
#Mean/median %: 3.54 1.64


if not np.isfinite(_scalar_objective_function(x0)):
    x0 = _unpack_theta(Theta())

op_params, code = op.leastsq(_objective_function, x0, maxfev=5000)

xp = _pack_theta(op_params)


Y_P = predict_line_strength(op_params)

rms = np.std(Y - Y_P)
percentage = 100 * np.abs((Y - Y_P)/Y)
print("RMS: {:.2f}".format(rms))
print("Mean/median %: {:.2f} {:.2f}".format(np.mean(percentage), np.median(percentage)))


for column_name in ("logg", "feh"):
    fig, ax = plt.subplots()
    scat = ax.scatter(Y, Y_P, c=data[column_name][only_consider], edgecolor="none")
    cbar = plt.colorbar(scat)
    cbar.set_label(column_name)


