import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

flow = []
flow.append(np.genfromtxt('../dump/flow1.txt'))
flow.append(np.genfromtxt('../dump/flow2.txt'))

# sns.set(context='notebook', style='darkgrid')
fig, ax = plt.subplots(1, 2, figsize=(9, 4))

x = np.linspace(0, 2, len(flow[0]))
# ax[0].plot(x, flow1[)
# ax[1].plot(x, flow2)
lbl = ['Азот', 'Пентан']

for k in 0, 1:
    flow[k] = flow[k].transpose()
    ax[k].set_title(lbl[k])
    # flow1 = (flow[k][0] + flow[k][1]) / 2
    flow2 = (flow[k][len(flow[k][0]) - 1]
           + flow[k][len(flow[k][0]) - 2]) / 2
    flow1 = flow[k][0]
    # flow2 = flow[k][len(flow[k][0]) - 1]
    ax[k].plot(x, flow1, label=r'$y = 0$')
    ax[k].plot(x, flow2, label=r'$y = 2$')
    ax[k].set_xlabel(r'$x, m$')
    ax[k].set_ylabel(r'$y, m$')
    ax[k].legend()

fig.savefig("../img/flows.pdf")
fig.savefig("../img/flows.png", dip=200)
plt.show()
