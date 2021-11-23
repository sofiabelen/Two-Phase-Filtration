import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


dt = [0.01, 0.005, 0.0001]
# dt = 1

sns.set(context='notebook', style='darkgrid')
fig, ax = plt.subplots(2, sharex=True, figsize=(8, 7))

for i in range(len(dt)):
    file = "../dump/mass" + str(dt[i]) + ".txt"
    data = np.genfromtxt(file, names=True)

    formula = r"$\sum_{i, j} \hat \rho_{i, j} \Delta x \Delta y + \sum_{i} \left( \Phi^{y}_{i, 1} - \Phi^{y}_{i, ny - 1} \right) \cdot dt \cdot \Delta x$"
    
    m1 = data['m1'] - data['m1'][0]
    m2 = data['m2'] - data['m2'][0]
    
    ax[0].plot(np.arange(0, len(m1)) * dt[i], m1, '--',
            label = "dt = " + str(dt[i]))
    ax[1].plot(np.arange(0, len(m2)) * dt[i], m2, '--',
            label = "dt = " + str(dt[i]))
    
ax[0].set_title("Nitrogen Mass: " + formula)
ax[1].set_title("Pentane Mass: "  + formula)
ax[1].set_xlabel(r'$t$')

for i in range(0, 2):
    ax[i].set_ylabel(r'$\Delta$mass')
    ax[i].legend()

fig.savefig("../img/mass-dt.png", dpi=200)
plt.show()
