import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
np.set_printoptions(precision=3)
from scipy.integrate import trapz
#import matplotlib
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'Times New Roman'
import matplotlib
matplotlib.font_manager._rebuild()
#mport matplotlib.font_manager
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'Times New Roman'
fm = matplotlib.font_manager.json_load("../../.cache/matplotlib/fontList.json")
fm.findfont("serif", rebuild_if_missing=False)
fm.findfont("serif", fontext="afm", rebuild_if_missing=False)

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

cmap = cm.RdPu_r
cmap2 = cm.viridis_r
sand = (0.9, 0.85, 0.)
mud = (0.66, 0.41, 0.21)
p_f = 0.5
sigma_T = 0.25


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def MandelRoots(N):
    eps = 0.00001
    zeros = np.zeros(N)
    for i in range(1,N):
        a1 = (i-1)*np.pi + np.pi/2 + eps
        a2 = a1 + np.pi - eps
        #print a1,a2
        for j in range(40):
            y1 = np.tan(a1)+a1/K
            y2 = np.tan(a2)+a2/K
            am = (a1+a2)/2.
            ym = np.tan(am)+am/K
            if (ym*y1>0): a1=am
            else: a2=am
            if (abs(y2)<eps): am=a2
        zeros[i] = am
        #print np.tan(zeros[i])+zeros[i]*K
    return zeros

def P_m(t):
    return t

def p_s(t,G):
    ps = np.zeros(len(t))
    p_f = 0.5
    dt = t[1]-t[0]
    for i in range(1,len(t)):
        ps[i] = ps[i-1] + G*dt
        if (ps[i] >= p_f):
            ps[i] = p_f - sigma_T
    return ps

def p_s_pipe(t,G):
    ps = np.zeros(len(t))
    p_f = 0.5
    dt = t[1]-t[0]
    for i in range(1,len(t)):
        ps[i] = ps[i-1] + G*dt
        if (ps[i] >= p_f):
            ps[i] = p_f - sigma_T
    return ps


def p_m_comp(x, t, N):
    zeros = MandelRoots(N)
    p = t + (G-1.) * ( (2*t - x*(2.-x))/(2.*(K+1.)) + K/(3.*(1. + K)**2) )
    for i in range(1,N):
        p -= 2.*K*(G-1.) * np.cos(zeros[i]*(1.-x)) / np.cos(zeros[i]) * np.exp(-zeros[i]**2*t) \
             / ( zeros[i]**2 * (K**2 + K + zeros[i]**2) )
    return p

def p_m_avg_comp(t, N):
    zeros = MandelRoots(N)
    p = t + (G-1.) * ( (2*t - 2/3.)/(2.*(K+1.)) + K/(3.*(1. + K)**2) )
    for i in range(1,N):
        p -= 2.*K*(G-1.) * np.tan(zeros[i]) * (np.exp(-zeros[i]**2*t)) \
             / ( zeros[i]**3 * (K**2 + K + zeros[i]**2) )
    return p

def p_s_comp(t, N):
    zeros = MandelRoots(N)
    p = t + (G-1.) * ( t/(K+1.) + K/(3.*(1. + K)**2) )
    for i in range(1,N):
        p -= 2.*K*(G-1.) * np.exp(-zeros[i]**2*t) \
             / ( zeros[i]**2 * (K**2 + K + zeros[i]**2) )
    return p

def p_m_pipe(x, t, N):
    zeros = MandelRoots(N)
    p = -sigma_T/(1. + K)
    for i in range(1,N):
        p -=  2.*K*sigma_T * np.cos(zeros[i]*(1.-x)) / np.cos(zeros[i]) * np.exp(-zeros[i]**2*t) \
             / ( K**2 + K + zeros[i]**2 )
    return p

def p_m_avg_pipe(t, t0, N):
    zeros = MandelRoots(N)
    p = -sigma_T/(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
            p -=  2.*K*sigma_T * np.tan(zeros[i]) * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( zeros[i] * ( K**2 + K + zeros[i]**2 ) ) * np.heaviside(t-t0, 1)
    return p

def p_s_pipe(t, t0, N):
    zeros = MandelRoots(N)
    p = -sigma_T/(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
        p -=  2.*K*sigma_T * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( K**2 + K + zeros[i]**2 ) * np.heaviside(t-t0, 1)
    return p

J = 101
x_grid = np.linspace(0, 1, J)
dx = x_grid[1]-x_grid[0]

lnw=1.25





G = 2.
K = 5.

Ts = 10000
t_grid = np.linspace(0, 2.8, Ts)
dt = t_grid[1]-t_grid[0]

sigma_T = 0.25
p_f = 0.5

alpha = dt/(2.*dx**2)
beta = K*dt/(2.*dx)

# initial condition
U = np.zeros(J)#+p_f
#U[0] = -1.

# matrix elements
A_u = np.diagflat([-alpha for i in range(J-1)], -1) +\
      np.diagflat([(1 + 2*alpha) + alpha/beta*(1 - beta)]+[1.+2.*alpha for i in range(J-2)]+[1.+alpha]) +\
      np.diagflat([-alpha]+[-alpha for i in range(J-2)], 1)

B_u = np.diagflat([alpha for i in range(J-1)], -1) +\
      np.diagflat([(1 - 2*alpha) + alpha/beta*(1 + beta)]+[1.-2.*alpha for i in range(J-2)]+[1.-alpha]) +\
      np.diagflat([alpha]+[alpha for i in range(J-2)], 1)

C_t = np.zeros(J)+1.
C_t[0] += G*alpha/beta
U_record = []
U_record.append(U)

avg_pm = np.zeros(len(t_grid))
ps = np.zeros(len(t_grid))

t_stop = 1.4

for ti in range(1,Ts):
    U_new = np.linalg.solve(A_u, B_u.dot(U)+dt*C_t)
    if (U_new[0] >= p_f):
        U_new[0] = p_f-sigma_T
        print ti*dt
    #if (ti*dt >= t_stop):
    #    C_t = np.zeros(J)
    U = U_new
    U_record.append(U)
    

for i in range(len(U_record)):
    avg_pm[i] = trapz(U_record[i], x_grid)
    ps[i] = U_record[i][0]

P_s = ps
P_m = avg_pm



fig, ax = plt.subplots(1, 1, figsize=(8,3))
idx = find_nearest(t_grid, t_stop)
ax.plot(t_grid[idx:], P_m[idx:], alpha=0.3, c='k', lw=1., ls=(0,(5,0.75)))











G = 2.
K = 5.

Ts = 10000
t_grid = np.linspace(0, 2.8, Ts)
dt = t_grid[1]-t_grid[0]

sigma_T = 0.25
p_f = 0.5

alpha = dt/(2.*dx**2)
beta = K*dt/(2.*dx)

# initial condition
U = np.zeros(J)#+p_f
#U[0] = -1.

# matrix elements
A_u = np.diagflat([-alpha for i in range(J-1)], -1) +\
      np.diagflat([(1 + 2*alpha) + alpha/beta*(1 - beta)]+[1.+2.*alpha for i in range(J-2)]+[1.+alpha]) +\
      np.diagflat([-alpha]+[-alpha for i in range(J-2)], 1)

B_u = np.diagflat([alpha for i in range(J-1)], -1) +\
      np.diagflat([(1 - 2*alpha) + alpha/beta*(1 + beta)]+[1.-2.*alpha for i in range(J-2)]+[1.-alpha]) +\
      np.diagflat([alpha]+[alpha for i in range(J-2)], 1)

C_t = np.zeros(J)+1.
C_t[0] += G*alpha/beta
U_record = []
U_record.append(U)

avg_pm = np.zeros(len(t_grid))
ps = np.zeros(len(t_grid))

t_stop = 1.4

for ti in range(1,Ts):
    U_new = np.linalg.solve(A_u, B_u.dot(U)+dt*C_t)
    if (U_new[0] >= p_f):
        U_new[0] = p_f-sigma_T
        print ti*dt
    if (ti*dt >= t_stop):
        C_t = np.zeros(J)
    U = U_new
    U_record.append(U)
    

for i in range(len(U_record)):
    avg_pm[i] = trapz(U_record[i], x_grid)
    ps[i] = U_record[i][0]

P_s = ps
P_m = avg_pm


ax.axhline(p_f, lw=0.75, c='k', alpha=1)
ax.axhline(p_f - sigma_T, lw=0.75, c='k', alpha=1)

ax.axvline(t_stop, lw=1.0, c='k', ls='dotted', alpha=1)

ax.plot(t_grid, G*t_grid, alpha=0.15, c='k', lw=1.)
ax.plot(t_grid, t_grid, alpha=0.3, c='k', lw=1., ls=(0,(5,0.75)))
ax.plot(t_grid, P_s, c=sand, lw=1., label='sandstone') #, ls=(0,(5,0.5)))
ax.plot(t_grid, P_m, c=mud, lw=1., label='mudstone', ls=(0,(5,0.75)))


ax.set_xlabel('time $t^*$')
ax.set_ylabel('average overpressure $\overline{p^*}$')

#xticks = np.arange(0, 1.7, 0.2)
#ax.set_xticks(xticks)

ax.set_xlim([0., t_grid[-1]])
ax.set_ylim([0., 1.])

ax.annotate('compression stops', xy=(0.51, 0.915), fontsize=12, textcoords='axes fraction', annotation_clip=False)

ax2 = ax.twinx()

ticks = [0.25, 0.5]
ticklabels = ['$p_f^* - \sigma_T^*$', '$p_f^*$']

ax2.set_yticks(ticks)
ax2.set_yticklabels(ticklabels)
ax2.set_frame_on(False)
plt.tight_layout()
plt.savefig('./sample.pdf', bbox_inches="tight")
#plt.show()

