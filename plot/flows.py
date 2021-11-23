import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

flow_in = []
flow_out = []
flow_in.append(np.genfromtxt('../dump/flow1_in.txt'))
flow_in.append(np.genfromtxt('../dump/flow2_in.txt'))
flow_out.append(np.genfromtxt('../dump/flow1_out.txt'))
flow_out.append(np.genfromtxt('../dump/flow2_out.txt'))

fig, ax = plt.subplots(1, 2, figsize=(9, 4))

x = np.linspace(0, 2, len(flow_in[0]))
lbl = ['Азот', 'Пентан']
formula = r'\rho \cdot v_y'

for k in 0, 1:
    # ax[k].plot(x, flow_in[k],  label=r'$y = 0$')
    ax[k].plot(x, flow_out[k], label=r'$y = 2$')
    ax[k].set_xlabel(r'$x, m$')
    ax[k].set_ylabel(r'$y, m$')
    ax[k].set_title(lbl[k])
    ax[k].legend()

fig.savefig("../img/flows.pdf")
fig.savefig("../img/flows.png", dip=200)
plt.show()
