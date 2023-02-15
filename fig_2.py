import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
cmap = cm.GnBu_r

# gamma
G = 2.
# nu (only affects the figure axes)
K = 5.

t = np.array([0.01, 0.25, 0.5, 0.75, 1.0]) 
fig, ax = plt.subplots(1,2, figsize=(6,2), gridspec_kw={'width_ratios': [2, 1]})

ax[0].axhline(0, color='gray', linewidth=1, alpha=0.5)
for j in range(len(t)):
	x = np.linspace(-1./K, 1., 10001)
	ps = np.zeros(len(x))
	for k in range(len(x)):
		if (x[k] >= 0.):
			ps[k] = t[j]
		else: ps[k] = G*t[j]
	ax[0].plot(ps, x, c=cmap(j/float(len(t))), lw=1.25)

ax[1].plot([0,1], [0,1], c=(0.66, 0.41, 0.21), lw=1.25, label='mudstone', ls=(0,(5,1)))
ax[1].plot([0,1], [0,2], c=(0.9, 0.85, 0.), lw=1.25, label='sandstone')

ax[0].set_xlim([0,G])
ax[0].set_ylim([-1./K,1])
ax[0].invert_yaxis()

ax[1].set_xlim([0,1])
ax[1].set_ylim([0,2])

ax[0].set_xlabel('overpressure $p^*$') 
ax[1].set_xlabel('time $t^*$')

ax[0].set_ylabel('depth $z$', labelpad=-10)
ax[1].set_ylabel('average overpressure $\overline{p^*}$')

ticks = [-1./K, 0, 1]
ticklabels = ["$-h_s$", "0", "$h_m$"]

ax[0].set_yticks(ticks)
ax[0].set_yticklabels(ticklabels)

plt.show()
