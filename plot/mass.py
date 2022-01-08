import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

data = np.genfromtxt('../dump/mass.txt', names=True)

dt = 0.01

sns.set(context='notebook', style='darkgrid')
fig, ax = plt.subplots(2, sharex=True, figsize=(8, 7))

formula = r"$\sum_{i, j} \hat \rho_{i, j} \Delta x \Delta y + \sum_{i} \left( \hat \rho_{i, 0} \cdot v_{i, 0} - \hat\rho_{i, ny} \cdot v_{i, ny}\right) \cdot dt \cdot \Delta x$"

m1 = data['m1'] # - data['m1'][0]
m2 = data['m2'] # - data['m2'][0]
phi1 = data['phi1']
phi2 = data['phi2']

ax[0].plot(np.arange(0, len(m1)) * dt, m1, '--',
        label = "dt = " + str(dt))
ax[1].plot(np.arange(0, len(m2)) * dt, m2, '--',
        label = "dt = " + str(dt))
# ax[0].plot(np.arange(0, len(phi1)) * dt, phi1, '--',
#         label = r"$\Phi_1$")
# ax[1].plot(np.arange(0, len(phi2)) * dt, phi2, '--',
#         label = r"$\Phi_2$")

ax[1].set_xlabel(r'$t$')
for i in range(0, 2):
    ax[i].set_ylabel(r'$\Delta$mass')
    ax[i].legend()

fig.savefig("../img/mass.png", dpi=200)
plt.show()
