import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

flow_in = []
flow_out = []
flow_in.append(np.genfromtxt('../dump-50/flow1_in.txt'))
flow_in.append(np.genfromtxt('../dump-50/flow2_in.txt'))
flow_out.append(np.genfromtxt('../dump-50/flow1_out.txt'))
flow_out.append(np.genfromtxt('../dump-50/flow2_out.txt'))

sns.set(context='talk', style='darkgrid')
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_box_aspect(1)

x = np.linspace(0, 2, len(flow_in[0]))
lbl = ['Азот', 'Пентан']
formula = r'\rho \cdot v_y'

for k in range(2):
    ax.plot(x, flow_out[k], label=lbl[k])

ax.set_ylim(-0.01, 0.165)
ax.set_xlabel(r'$x, m$')
ax.set_ylabel(r'$\rho \cdot v_y, \frac{kg}{m^2 \cdot s}$')
ax.legend()

fig.savefig("../img/flows.pdf")
fig.savefig("../img/flows.png", dpi=200)
plt.show()
