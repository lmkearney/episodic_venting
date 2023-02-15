import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import trapz
import matplotlib
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

cmap2 = cm.viridis_r

def roots(K, N):
    eps = 0.00001
    zeros = np.zeros(N)
    for i in range(1,N):
        a1 = (i-1)*np.pi + np.pi/2 + eps
        a2 = a1 + np.pi - eps
        for j in range(40):
            y1 = K*np.tan(a1)+a1
            y2 = K*np.tan(a2)+a2
            am = (a1+a2)/2.
            ym = K*np.tan(am)+am
            if (ym*y1>0): a1=am
            else: a2=am
            if (abs(y2)<eps): am=a2
        zeros[i] = am
    return zeros

def periodic_solution(t, tau, K):
    zeros = roots(K, N)
    p = t/(1+K)
    for k in range(100):
	for i in range(1,N):
            p += 2. * K * (1.-np.exp(-zeros[i]**2*t * tau)) \
                / ( K**2 + K + zeros[i]**2 ) * np.exp(-zeros[i]**2*k*tau)
    return p 


fig, ax = plt.subplots(1,2, figsize=(6,3))
time = np.linspace(0, 1, 10001)

p_f = 0.
sigma_T = 0.1

# varying period
N = 100
K = 5.

period = 0.01
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0))

period = 0.1
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0.25/2.))

period = 0.5
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0.5/2.))

N = 1000 # solutions with large period require more terms for accuracy
period = 2.5
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0.75/2.))

ax[0].set_title('$\\nu = 5$')
ax[0].set_xlim([0,1])
ax[0].set_ylim([0,1])
ax[0].set_xlabel('time $t^*/\\tau^*$')
ax[0].set_ylabel('sandstone overpressure $p_s^*/\sigma_T^*$')
ax[0].set_aspect(1.0/ax[0].get_data_ratio(), adjustable='box')

# varying nu
N = 100
period = 5. 

K = 0.2
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0), label='$\\nu = 0.2$')

K = 1.
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0.25/2.), label='$\\nu = 1$')

K = 5.
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0.5/2.), label='$\\nu = 5$')

N = 2000 # solutions with large nu require more terms for accuracy
K = 50.
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0.75/2.), label='$\\nu = 50$')

ax[0].plot([0,1], [1.-1./(1.+5.), 1.], c='k', lw=1.)
ax[0].plot([0,1], [0,1], c='k', linestyle='--', lw=1.)
ax[1].plot([0,1], [0,1], c='k', linestyle='--', lw=1.)

ax[1].set_xlabel('time $t^*/\\tau^*$')
ax[1].set_title('$\\tau^* = 5$')
ax[1].set_xlim([0,1])
ax[1].set_ylim([0,1])

ax2 = ax[1].twinx()
ticks = [0,1]
ticklabels = ['$\displaystyle{\\frac{p_f^*}{\sigma_T^*} - 1}$', '$\displaystyle{\\frac{p_f^*}{\sigma_T^*}}$']
ax2.set_ylim([0,1])
ax2.set_yticks(ticks)
ax2.set_yticklabels(ticklabels)
ax2.set_frame_on(False)

plt.tight_layout()
plt.show()

