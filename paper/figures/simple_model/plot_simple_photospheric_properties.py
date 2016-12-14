
"""
Plot photospheric properties as a function of optical depth for Solar-metallicity
starts.
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, rc
from matplotlib.ticker import MaxNLocator
from smh.photospheres import interpolator as PhotosphereInterpolator



isochrone_filename = "dartmouth-solar-isochrone.iso"

logTeffs, loggs = np.loadtxt(isochrone_filename, usecols=(2, 3), unpack=True)
teffs = 10**logTeffs

photosphere_interpolator = PhotosphereInterpolator()


# Interpolate stellar photospheres.
fig, axes = plt.subplots(2)

for teff, logg in zip(teffs, loggs):
    if not 6500 > teff > 4000: continue

    photosphere = photosphere_interpolator(teff, logg, 0.0)

    x = np.log10(photosphere["RHOX"])
    y1 = photosphere["T"]
    y2 = np.log10(photosphere["ABROSS"])

    axes[0].plot(x, y1, alpha=0.01, c="#666666")
    axes[1].plot(x, y2, alpha=0.01, c="#666666")


# Set common x-axes.
limits = np.array([ax.get_xlim() for ax in axes])
limits = (np.min(limits), np.max(limits))
for ax in axes:
    ax.set_xlim(limits)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.set_aspect(np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()), adjustable="box")

axes[0].set_xticklabels([])
axes[1].set_xlabel(r"Optical depth, $\log_{10}\tau$")

axes[0].set_ylabel(r"Temperature, $T$ $({\rm K})$")
axes[1].set_ylabel(r"Rosseland opacity, $\kappa$")

fig.tight_layout()
fig.subplots_adjust(left=0.15, bottom=0.15)


fig.savefig("simple_photospheric_properties.pdf", dpi=300, bbox_inches="tight")



