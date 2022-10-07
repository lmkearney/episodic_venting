import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
matplotlib.font_manager._rebuild()
#mport matplotlib.font_manager
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'Times New Roman'
fm = matplotlib.font_manager.json_load("../../.cache/matplotlib/fontList.json")
fm.findfont("serif", rebuild_if_missing=False)
fm.findfont("serif", fontext="afm", rebuild_if_missing=False)
sand = (0.9, 0.85, 0.)
mud = (0.66, 0.41, 0.21)

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
    ps = p_0 + np.zeros(len(t))
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

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

G_nu = 1. + 10**np.linspace(-3, 3, 10000) #np.linspace(0, 1000, 10000) 
Contr = 1./(G_nu)

fig, ax = plt.subplots(figsize=(5.55,3))


ax.axvline(2., c='k', alpha=0.5, lw=0.75, ls='dotted')


ax.semilogx(G_nu, Contr, c='orangered', label='compression')
ax.set_xlim([1e0, 1e2])
ax.set_ylim([0,1])
#ax[0].set_ylabel('period $\\tau \, / \, \\tau_\mathrm{comp}$')
ax.set_ylabel('period $\displaystyle \\frac{\\tau \Gamma_s}{\sigma_T}$')
#ax[0].set_xlabel('$\\nu/\gamma = h_m/h_s \, \\alpha_m/\\alpha_s \, v_s/v_m$')
#ax[0].set_xlabel('$h_m \\alpha_m (1 + \lambda_s)/ h_s \\alpha_s (1 + \lambda_m)$')
ax.set_xlabel('venting frequency multiplier $\displaystyle 1 + \\frac{\\nu}{\gamma} = 1 + \\frac{h_m}{h_s} \, \\frac{\\alpha_m}{\\alpha_s} \, \\frac{1 + \lambda_s/\mu_s}{1 + \lambda_m/\mu_m}$')
#plt.legend()

#ax[1].set_aspect(1.0/ax[1].get_data_ratio(), adjustable='box')
ax2 = ax.twinx()

ticks = np.linspace(0,1,6)
ax2.set_yticks(ticks)
ax2.set_ylabel('pressure recharge contribution \n from compression')
#ax.set_aspect('equal')#1.0/ax.get_data_ratio(), adjustable='box')

ax.annotate('compression-\ndominated', xy=(-0.02, 1.05), fontsize=12, textcoords='axes fraction', annotation_clip=False)
ax.annotate('diffusion-dominated', xy=(0.4, 1.075), fontsize=12, textcoords='axes fraction', annotation_clip=False)
#ax.annotate('$\displaystyle{\\frac{\\nu}{\gamma} = 1}$', xy=(0.15,-0.15), fontsize=10, textcoords='axes fraction', annotation_clip=False)
ax.annotate('$\displaystyle{\\frac{\\nu}{\gamma} = 1}$', xy=(0.165,0.04), fontsize=10, textcoords='axes fraction', annotation_clip=False)


#ax.annotate('compression \n dominated', xy=(0.09, 0.325), fontsize=12, textcoords='axes fraction', annotation_clip=False)
#ax.annotate('diffusion \n dominated', xy=(0.7, 0.545), fontsize=12, textcoords='axes fraction', annotation_clip=False)
J = 101
x_grid = np.linspace(0, 1, J)
dx = x_grid[1]-x_grid[0]

N = 501
t_grid = np.linspace(0, 2.0, N)
dt = t_grid[1]-t_grid[0]

lnw=1.25
G = 2.
K = 8.
time = np.linspace(0, 2.0, 10001) #10**t
N=1000

sigma_T = 1.
p_0 = 0.9
p_f = 1.
axins = inset_axes(ax, width=1, height=0.7, loc="lower left", 
		   bbox_to_anchor=(0.2,0.55,0,0), bbox_transform=ax.transAxes)
axins2 = inset_axes(ax, width=1, height=0.7, loc="lower left",
		    bbox_to_anchor=(0.6,0.15,0,0), bbox_transform=ax.transAxes)

axins.set_xticks([])
axins.set_yticks([])
axins.set_xlabel('time')
axins.set_ylabel('$p_s$')
axins.set_xlim([1,2])
axins.set_ylim([0,1])
#axins.xaxis.set_label_position('top') 
axins.plot(time, p_s(time, G), c=sand)

axins2.set_xticks([])
axins2.set_yticks([])
axins2.set_xlabel('time')
axins2.set_ylabel('$p_s$')

P_s = p_0+p_s_comp(time, N)
P_m = p_0+p_m_avg_comp(time, N)

t_event = 0
c = 0
while (any(P_s>p_f)):
    t_event = time[np.argmax(P_s>p_f)]
    print t_event
    P_s += p_s_pipe(time, t_event, N)
    P_m += p_m_avg_pipe(time, t_event, N)
    c+=1
axins2.plot(time, P_s, c=sand)
axins2.set_xlim([1,2])
axins2.set_ylim([0,1])
print sigma_T/G/(1.+K/G)
plt.tight_layout()

#ax.set_zorder(0)
#ax2.set_zorder(1)
ax.plot(1. + K/G, 1./(1.+K/G), 'o', zorder=10, clip_on=False, c='royalblue')
ax.plot(1,1, 'o', zorder=10, clip_on=False, c='seagreen')

axins.plot(1, 0, 'o', zorder=10, clip_on=False, c='seagreen', ms=4.)
axins2.plot(1, 0, 'o', zorder=10, clip_on=False, c='royalblue', ms=4.)
#axins2.tick_params(color='royalblue', labelcolor='k')
#for spine in axins2.spines.values():
#    spine.set_edgecolor('royalblue')

#axins2.spines['top'].set_visible(False)
#axins2.spines['right'].set_visible(False)
#axins2.spines['bottom'].set_visible(False)
#axins2.spines['left'].set_visible(False)
#ax.minorticks_off()
ax2.set_frame_on(False)

plt.savefig('./t_vs_d2.pdf', bbox_inches="tight")
plt.show()


