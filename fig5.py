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

def p_m_pipe(x, t, G, K):
    zeros = roots(K, N)
    p = -1./(1. + K)
    for i in range(1,N):
        p -=  2.*K * np.cos(zeros[i]*(1.-x)) / np.cos(zeros[i]) * np.exp(-zeros[i]**2*t) \
             / ( K**2 + K + zeros[i]**2 )
    return p

def p_m_avg_pipe(t, t0, G, K):
    zeros = roots(K, N)
    p = -1./(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
            p -=  2.*K * np.tan(zeros[i]) * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( zeros[i] * ( K**2 + K + zeros[i]**2 ) ) * np.heaviside(t-t0, 1)
    return p

def p_s_pipe(t, t0, G, K):
    zeros = roots(K, N)
    p = -1./(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
        p -=  2.*K * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( K**2 + K + zeros[i]**2 ) * np.heaviside(t-t0, 1)
    return p

N = 500
J = 1001
x = np.linspace(0, 1, J)

Ks = [0.2, 1., 5.]
G = 2.

t = np.array([0.01, 0.25, 0.50, 0.75, 1.00])
fig, ax = plt.subplots(3,2, figsize=(5,6), gridspec_kw={'width_ratios': [2, 1]})

for i in range(3):
    ax[i][0].axhline(0, color='gray', linewidth=1, alpha=0.5)
    for j in range(len(t)):
        ax[i][0].plot(p_m_pipe(x, t[j], G, Ks[i]), x, c=cmap(j/(float(len(t)))), lw=1.25)
	x_2 = np.linspace(0, -1./Ks[i], J)
        ax[i][0].plot(np.zeros(len(x_2)) + p_m_pipe(0, t[j], G, Ks[i]), x_2, c=cmap(j/(float(len(t)))), lw=1.25)
    ax[i][0].axvline(-1/(1+Ks[i]), c='k', lw=0.75, alpha=0.5)
    
    ax[i][1].axhline(-1/(1+Ks[i]), c='k', lw=0.75, alpha=0.5)
    ax[i][1].plot(x, p_s_pipe(x, 0., G, Ks[i]), c=(0.9, 0.85, 0.), lw=1.25, label='sandstone')
    ax[i][1].plot(x, p_m_avg_pipe(x, 0., G, Ks[i]), c=(0.66, 0.41, 0.21), lw=1.25, ls=(0,(5,1)), label='mudstone')


ticklabels = ["$-h_s$", "0", "$h_m$"]
for i in range(3):
    ticks = [-1./Ks[i], 0, 1]
    ax[i][0].set_yticks(ticks)
    ax[i][0].set_yticklabels(ticklabels)
    ax[i][0].set_xlim([-1., 0])
    ax[i][0].set_ylim([-1./Ks[i],1])
    ax[i][0].invert_yaxis()
    ax[i][0].set_aspect(1./ax[i][0].get_data_ratio(), adjustable='box')

    ax[i][1].set_xlim([0,1])
    ax[i][1].set_ylim([-1,0])
    ax[i][1].set_aspect(1./ax[i][1].get_data_ratio(), adjustable='box')

ax[0][0].set_ylabel('depth $z$', labelpad=-10)
ax[0][1].set_ylabel('average overpressure $\overline{p}/\sigma_T$')

ax[2][1].set_xlabel('time $t^*$')
ax[2][0].set_xlabel('overpressure $p/\sigma_T$')

plt.tight_layout()
plt.show()

