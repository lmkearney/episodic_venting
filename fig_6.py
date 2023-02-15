import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

cmap = cm.RdPu_r
sand = (0.9, 0.85, 0.)
mud = (0.66, 0.41, 0.21)

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

def p_s_sawtooth(t,G):
    ps = np.zeros(len(t))
    p_f = 0.5
    dt = t[1]-t[0]
    for i in range(1,len(t)):
        ps[i] = ps[i-1] + G*dt
        if (ps[i] >= p_f):
            ps[i] = p_f - sigma_T
    return ps

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

def p_s_pipe(t,G):
    ps = np.zeros(len(t))
    p_f = 0.5
    dt = t[1]-t[0]
    for i in range(1,len(t)):
        ps[i] = ps[i-1] + G*dt
        if (ps[i] >= p_f):
            ps[i] = p_f - sigma_T
    return ps

def p_m_avg_pipe(t, t0, K):
    zeros = roots(K, N)
    p = -sigma_T/(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
            p -=  2.*K*sigma_T * np.tan(zeros[i]) * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( zeros[i] * ( K**2 + K + zeros[i]**2 ) ) * np.heaviside(t-t0, 1)
    return p

def p_s_pipe(t, t0, K):
    zeros = roots(K, N)
    p = -sigma_T/(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
        p -=  2.*K*sigma_T * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( K**2 + K + zeros[i]**2 ) * np.heaviside(t-t0, 1)
    return p

N = 1000 # number of terms in infinite series

p_f = 0.5
sigma_T = 0.25
G = 2. # gamma
K = 5. # nu

J = 10001
time = np.linspace(0, 2.0, J)

# initial solution without venting while P_s < p_f
P_s = p_s_comp(time, G, K)
P_m = p_m_avg_comp(time, G, K)

# add venting solutions whenever P_s > p_f
t_event = 0
c = 0
while (any(P_s>p_f)):
    t_event = time[np.argmax(P_s>p_f)]
    P_s += p_s_pipe(time, t_event, K)
    P_m += p_m_avg_pipe(time, t_event, K)
    c+=1

fig, ax = plt.subplots(1, 2, figsize=(8,3))

ax[0].axhline(p_f, lw=0.75, c='k', alpha=1)
ax[0].axhline(p_f - sigma_T, lw=0.75, c='k', alpha=1)
ax[1].axhline(p_f, lw=0.75, c='k', alpha=1)
ax[1].axhline(p_f - sigma_T, lw=0.75, c='k', alpha=1)

# compression only plots
ax[0].plot(time, G*time, alpha=0.15, c='k', lw=1.)
ax[0].plot(time, p_s_sawtooth(time, G), c=sand, lw=1., label='sandstone')
ax[0].plot(time, time, c=mud, lw=1., ls=(0,(5,0.75)), label='mudstone')

# compression and diffusion plots
ax[1].plot(time, G*time, alpha=0.15, c='k', lw=1.)
ax[1].plot(time, time, alpha=0.3, c='k', lw=1., ls=(0,(5,0.75)))
ax[1].plot(time, P_s, c=sand, lw=1., label='sandstone') 
ax[1].plot(time, time, c=mud, lw=1., label='mudstone', ls=(0,(5,0.75)))

ax[0].set_xlabel('time $t^*$')
ax[1].set_xlabel('time $t^*$')
ax[0].set_ylabel('average overpressure $\overline{p^*}$')

xticks = np.arange(0, 1.7, 0.2)
ax[0].set_xticks(xticks)
ax[1].set_xticks(xticks)

ax[0].set_xlim([0., 1.4])
ax[1].set_xlim([0., 1.4])

ax[0].set_ylim([0., 1.])
ax[1].set_ylim([0., 1.])

ax[0].set_title('compression only')
ax[1].set_title('compression and diffusion')
ax2 = ax[1].twinx()

# labelling p_f^* and p_f^*-\sigma_T^* on the axes
ticks = [0.25, 0.5]
ticklabels = ['$p_f^* - \sigma_T^*$', '$p_f^*$']
ax2.set_yticks(ticks)
ax2.set_yticklabels(ticklabels)
ax2.set_frame_on(False)

plt.tight_layout()
plt.show()
