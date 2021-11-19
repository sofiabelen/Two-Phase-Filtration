import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

data = np.genfromtxt('../dump/mass.txt', names=True)

dt = 0.0001
# dt = 1

sns.set(context='notebook', style='darkgrid')
fig, ax = plt.subplots(2, sharex=True, figsize=(8, 7))

formula = r"$\sum_{i, j} \hat \rho_{i, j} \Delta x \Delta y + \sum_{i} \left( \hat \rho_{i, 0} \cdot v_{i, 0} - \hat\rho_{i, ny} \cdot v_{i, ny}\right) \cdot dt \cdot \Delta x$"

ax[0].plot(np.arange(0, len(data['m1'])) * dt, data['m1'], '--',
        label = "Nitrogen Mass: " + formula + "\n dt = " + str(dt))
ax[1].plot(np.arange(0, len(data['m2'])) * dt, data['m2'], '--',
        label = "Pentane Mass: " + formula  + "\n dt = " + str(dt))

ax[1].set_xlabel(r'$t$')
for i in range(0, 2):
    ax[i].set_ylabel(r'$mass$')
    ax[i].legend()

fig.savefig("../img/mass.png", dpi=200)
plt.show()
