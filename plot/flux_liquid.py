import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

dt = 0.01
influx  = np.genfromtxt('../dump/total_flux_init.txt',  names=True)
outflux = np.genfromtxt('../dump/total_flux_final.txt', names=True)

fig, ax = plt.subplots(figsize=(7, 7))

t = np.arange(0, len(influx['phi1'])) * dt

ax.plot(t, influx['phi1'],  label=r'$y = 0$')
ax.plot(t, outflux['phi1'], label=r'$y = 2$')
ax.set_xlabel(r'$t, s$')
ax.set_ylabel(r'$mass, kg$')
ax.legend()

fig.savefig("../img/total_flux_liquid.pdf")
fig.savefig("../img/total_flux_liquid.png", dpi=200)
plt.show()
