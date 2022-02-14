import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

dt = 0.01
steps = 100000

t = 1000

flow_in = []
flow_out = []
flow_in.append(np.genfromtxt('../dump-' + str(t) + '/flow1_in.txt'))
flow_in.append(np.genfromtxt('../dump-' + str(t) + '/flow2_in.txt'))
flow_out.append(np.genfromtxt('../dump-' + str(t) + '/flow1_out.txt'))
flow_out.append(np.genfromtxt('../dump-' + str(t) + '/flow2_out.txt'))

sns.set(context='notebook', style='darkgrid')
fig, ax = plt.subplots(1, 2, figsize=(9, 5))

x = np.linspace(0, 2, len(flow_in[0]))
lbl = ['Азот', 'Пентан']
# lbl = ['Nitrogen', 'Pentane']
formula = r'\rho \cdot v_y'
c = ['', '']

for k in 0, 1:
    # ax[k].plot(x, flow_in[k],  label=r'$y = 0$')
    ax[k].plot(x, flow_out[k], label=r'$y = 2$')
    ax[k].set_xlabel(r'$x, m$')
    ax[k].set_title(lbl[k])
    ax[k].legend()

ax[0].set_ylim(-0.01, 0.165)
ax[1].set_ylim(-0.001, 0.0065)

# fig.suptitle(r"Потоки через выход в момент времени $t = $"
#         + str(int(steps * dt)))
# fig.suptitle(r"Fluxes through the exit at time $t = $"
#         + str(int(t)))
ax[0].set_ylabel(r'$\rho \cdot v_y, \frac{kg}{m \cdot s}$')
# fig.savefig("../../img/flows.pdf")
# fig.savefig("../../img/flows.png", dpi=200)
# fig.savefig("../../sofiabelen.github.io/images/flows-" + str(t) + ".png", dpi=200)
fig.savefig("../../Two-Phase-Filtration/gos-exam/img/flows-" + str(t) + ".pdf")
# plt.show()
