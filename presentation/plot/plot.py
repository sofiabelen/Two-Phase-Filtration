import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set(context='talk', style='whitegrid')
fig, ax = plt.subplots(figsize=(7, 7))

ax.set_xlabel(r'$s$')
ax.set_ylabel(r'$f(s)$')

x = np.arange(0, 1, 0.01)
ax.plot(x, x**2, '--',       label=r'$f_1(s)$ - удельная проницаемость газа')
ax.plot(x, (1 - x)**2, '--', label=r'$f_2(s)$ - удельная проницаемость жидкости')

ax.legend()

fig.savefig("../img/two-phase.pdf")
plt.show()
