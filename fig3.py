import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

cmap = cm.GnBu_r

def roots(K, N):
    eps = 0.00001
    zeros = np.zeros(N)
    for i in range(1,N):
        a1 = (i-1)*np.pi + np.pi/2 + eps
        a2 = a1 + np.pi - eps
        
        for j in range(40):
            y1 = np.tan(a1)+a1/K
            y2 = np.tan(a2)+a2/K
            am = (a1+a2)/2.
            ym = np.tan(am)+am/K
            if (ym*y1>0): a1=am
            else: a2=am
            if (abs(y2)<eps): am=a2
        zeros[i] = am
        
    return zeros

def p_m_comp(x, t, G, K):
    zeros = roots(K, N)
    p = t + (G-1.) * ( (2*t - x*(2.-x))/(2.*(K+1.)) + K/(3.*(1. + K)**2) )
    for i in range(1,N):
        p -= 2.*K*(G-1.) * np.cos(zeros[i]*(1.-x)) / np.cos(zeros[i]) * np.exp(-zeros[i]**2*t) \
             / ( zeros[i]**2 * (K**2 + K + zeros[i]**2) )
    return p

def p_m_avg_comp(t, G, K):
    zeros = roots(K, N)
    p = t + (G-1.) * ( (2*t - 2/3.)/(2.*(K+1.)) + K/(3.*(1. + K)**2) )
    for i in range(1,N):
        p -= 2.*K*(G-1.) * np.tan(zeros[i]) * (np.exp(-zeros[i]**2*t)) \
             / ( zeros[i]**3 * (K**2 + K + zeros[i]**2) )
    return p

def p_s_comp(t, G, K):
    zeros = roots(K, N)
    p = t + (G-1.) * ( t/(K+1.) + K/(3.*(1. + K)**2) )
    for i in range(1,N):
        p -= 2.*K*(G-1.) * np.exp(-zeros[i]**2*t) \
             / ( zeros[i]**2 * (K**2 + K + zeros[i]**2) )
    return p

N = 1000
J = 501
x_grid = np.linspace(0, 1, J)
dx = x_grid[1]-x_grid[0]

t_grid = np.linspace(0, 1., J)
dt = t_grid[1]-t_grid[0]

Ks = [0.2, 1., 5.]
G = 2.
t = np.linspace(0.01,1,5)
fig, ax = plt.subplots(3,2, figsize=(6,6), gridspec_kw={'width_ratios': [2, 1]})

for i in range(3):
    ax[i][1].plot([0,1], [0,1], c='k', lw=1.5, alpha=0.3, ls=(0,(5,1)))
    ax[i][1].plot([0,1], [0,2], c='k', lw=1.5, alpha=0.15)
    ax[i][1].set_xlim([0,1])
    ax[i][1].set_ylim([0,2])

    ax[i][0].axhline(0, color='gray', linewidth=1, alpha=0.5)
    for j in range(len(t)):
        ax[i][0].plot(p_m_comp(x_grid, t[j], G, Ks[i]), x_grid, c=cmap(j/(float(len(t)))), lw=1.5)
	x_2 = np.linspace(0, -1./Ks[i], J)
        ax[i][0].plot(np.zeros(len(x_2)) + p_m_comp(0, t[j], G, Ks[i]), x_2, c=cmap(j/(float(len(t)))), lw=1.5)

    ax[i][1].plot(t_grid, p_s_comp(t_grid, G, Ks[i]), c=(0.9, 0.85, 0.), lw=1.5, label='sandstone')
    ax[i][1].plot(t_grid, p_m_avg_comp(t_grid, G, Ks[i]), c=(0.66, 0.41, 0.21), lw=1.5, label='mudstone', ls=(0,(5,1)))

ticklabels = ["$-h_s$", "0", "$h_m$"]
titles = ['0.2', '1', '5']
for i in range(3):
    ticks = [-1./Ks[i], 0, 1]
    ax[i][0].set_yticks(ticks)
    ax[i][0].set_yticklabels(ticklabels)
    ax[i][0].set_xlim([0,G])
    ax[i][0].set_ylim([-1./Ks[i],1])
    ax[i][0].invert_yaxis()

ax[0][0].set_ylabel('depth $z$', labelpad=-10)
ax[0][1].set_ylabel('average overpressure $\overline{p^*}$')

ax[2][1].set_xlabel('time $t^*$')
ax[2][0].set_xlabel('overpressure $p^*$')

ax[0][0].annotate('$\\nu = 0.2$', xy=(-0.4, 0.45), fontsize=14, textcoords='axes fraction', annotation_clip=False)
ax[1][0].annotate('$\\nu = 1$', xy=(-0.4, 0.45), fontsize=14, textcoords='axes fraction', annotation_clip=False)
ax[2][0].annotate('$\\nu = 5$', xy=(-0.4, 0.45), fontsize=14, textcoords='axes fraction', annotation_clip=False)

plt.tight_layout()
plt.show()
