
"""
Plot the effective temperatures and surface gravities for Solar-metallicity
Dartmouth isochrones.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, rc


isochrone_filename = "dartmouth-solar-isochrone.iso"

logTeff, logg = np.loadtxt(isochrone_filename, usecols=(2, 3), unpack=True)
teff = 10**logTeff

# Parse ages.
with open(isochrone_filename, "r") as fp:
    contents = fp.readlines()

age = np.nan
ages = []
for line in contents:
    if line.startswith("#AGE"):
        age = float(line.split("=")[1].lstrip().split(" ")[0])
        continue

    if line.count(" ") > 13 and not line.startswith("#EEP") and np.isfinite(age):
        ages.append(age)

ages = np.array(ages)
assert ages.size == teff.size


fig, ax = plt.subplots()

cmap = cm.get_cmap("bone_r", len(set(ages)))


styles = {
    1.0: ("#93D0BF", 100),
    2.5: ("#72BFC4", 75), 
    5.0: ("#52AEC9", 50),
    7.5: ("#4095B5", 30),
    10.0: ("#3B748A", 20),
    12.0: ("#37535E", 10)
}
for age in sorted(list(set(ages))):
    match = ages == age

    #f = (age + 5)/20.0
    color, size = styles[age]
    ax.scatter(teff[match], logg[match],
        s=size,
        facecolor=color, edgecolor=color, linewidth=0.5, alpha=0.50,
        label=r"{:.1f} Gyr".format(age))

ax.set_xlim(8500, 3000)
ax.set_ylim(5.5, -0.5)

legend = plt.legend(loc="upper left", frameon=False, fontsize=12, borderpad=1.0,
    numpoints=1, scatterpoints=1)

"""
renderer = fig.canvas.get_renderer()
shift = max([t.get_window_extent(renderer).width for t in legend.get_texts()])
for t in legend.get_texts():
    t.set_ha('right') # ha is alias for horizontalalignment
    t.set_position((shift,0))
"""

ax.set_xlabel(r"Effective temperature, $T_{\rm eff}$ $({\rm K})$")
ax.set_ylabel(r"Surface gravity, $\log{g}$")

ax.set_aspect(np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim()), adjustable="box")

fig.tight_layout()
fig.subplots_adjust(left=0.15, bottom=0.15)

fig.savefig("simple_solar_isochrone.pdf", dpi=300, bbox_inches="tight")



