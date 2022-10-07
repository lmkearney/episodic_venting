import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import trapz
import matplotlib
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

cmap = cm.RdPu_r

p_f = 0.
sigma_T = 0.1

def MandelRoots(N):
    eps = 0.00001
    zeros = np.zeros(N)
    for i in range(1,N):
        a1 = (i-1)*np.pi + np.pi/2 + eps
        a2 = a1 + np.pi - eps
        #print a1,a2
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

def p_m_avg_pipe(t, t0, G, K):
    zeros = MandelRoots(N)
    p = -sigma_T/(1. + K) * np.heaviside(t-t0, 1)
    for i in range(1,N):
            p -=  2.*K*sigma_T * np.tan(zeros[i]) * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( zeros[i] * ( K**2 + K + zeros[i]**2 ) ) * np.heaviside(t-t0, 1)
    return p

def p_s_pipe(t, t0, K):
    zeros = MandelRoots(N)
    p = -sigma_T/(1. + K)
    for i in range(1,N):
        p_i = 2.*K*sigma_T * np.exp(-zeros[i]**2*abs(t-t0)) \
                / ( K**2 + K + zeros[i]**2 )
        p -= p_i
        #print i, p_i
    return p

def periodic_solution(t, tau, K):
    zeros = MandelRoots(N)
    p = t/(1+K)
    for k in range(100):
	for i in range(1,N):
            p += 2. * K * (1.-np.exp(-zeros[i]**2*t * tau)) \
                / ( K**2 + K + zeros[i]**2 ) * np.exp(-zeros[i]**2*k*tau)
    return p 

N = 100

K = 5.
time = np.linspace(0, 1, 10001)

fig, ax = plt.subplots(1,2, figsize=(6,3))

period = 0.01 #sigma_T/(K+G)
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0), label='$\gamma = 0.2$')

period = 0.1 #sigma_T/(K+G)
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0.25/2.), label='$\gamma = 1$')

period = 0.5 #sigma_T/(K+G)
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0.5/2.), label='$\gamma = 5$')

N = 1000
period = 2.5
ax[0].plot(time, periodic_solution(time, period, K), c=cmap2(0.75/2.), label='$\gamma = 50$')

#period = 0.
#ax[0].plot(time, periodic_solution(time, period, K), c='k', linestyle=(0,(5,1)), lw=1.)

#ax[0].plot([0,1], [0,1], 'k--', lw=1.)

ax[0].set_title('$\\nu = 5$')
ax[0].set_xlim([0,1])
ax[0].set_ylim([0,1])
ax[0].set_xlabel('time $t^*/\\tau^*$')
ax[0].set_ylabel('sandstone overpressure $p_s^*/\sigma_T^*$')
#ax[0].legend()

ax[0].set_aspect(1.0/ax[0].get_data_ratio(), adjustable='box')



N = 100
K = 0.2
period = 5. 
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0), label='$\\nu = 0.2$')

K = 1.
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0.25/2.), label='$\\nu = 1$')

K = 5.
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0.5/2.), label='$\\nu = 5$')

N = 2000
K = 50.
ax[1].plot(time, periodic_solution(time, period, K), c=cmap2(0.75/2.), label='$\\nu = 50$')

ax[0].plot([0,1], [1.-1./(1.+5.), 1.], c='k', lw=1.)
ax[0].plot([0,1], [0,1], c='k', linestyle='--', lw=1.)
ax[1].plot([0,1], [0,1], c='k', linestyle='--', lw=1.)
ax[1].set_xlabel('time $t^*/\\tau^*$')

#leg = ax[0].legend(loc='lower right')
#leg.get_frame().set_facecolor('none')
#leg.get_frame().set_linewidth(0.0)

#leg = ax[1].legend(loc='lower right')
#leg.get_frame().set_facecolor('none')
#leg.get_frame().set_linewidth(0.0)
#ax[1].legend()
ax[1].set_title('$\\tau^* = 5$')
ax[1].set_xlim([0,1])
ax[1].set_ylim([0,1])
#ax[1].set_xlabel('time $t/\\tau$')
#ax[1].set_ylabel('sandstone pressure $p_s/\sigma_T$')
#ax[0].arrow(0.2, 0.62, 0.18, -0.18, overhang=0.4, head_width=0.035, head_length=0.05, lw=1., ec='k', fc='k', zorder=100, clip_on=False)
#ax[1].arrow(0.4, 0.4, -0.3, 0.4, overhang=0.4, head_width=0.035, head_length=0.05, lw=1., ec='k', fc='k', zorder=100, clip_on=False)


#import matplotlib.patches as patches
#style = "fancy, tail_width=0.1, head_width=7.5, head_length=15"
#kw = dict(arrowstyle=style, color="k")
#
#a2 = patches.FancyArrowPatch((0.9, 0.8), (0.8, 1.1),
#                             connectionstyle="arc3,rad=-.4", zorder=100, clip_on=False, lw=1, mutation_scale=0.4, capstyle="butt", **kw)
#
#a3 = patches.FancyArrowPatch((0.9, 0.8), (0.8, 1.1),
#                             connectionstyle="arc3,rad=-.4", zorder=100, clip_on=False, lw=1, mutation_scale=0.4, **kw)
#
#a3.set_joinstyle('miter')
#
#ax[0].add_patch(a2)
#ax[1].add_patch(a3)

#ax[0].annotate('$\\tau^*$', xy=(0.825, 1.07), fontsize=14, textcoords='axes fraction', annotation_clip=False)
#ax[1].annotate('$\\nu$', xy=(0.825, 1.07), fontsize=14, textcoords='axes fraction', annotation_clip=False)


ax[0].annotate('$\\tau^* \\gg 1$', xy=(0.05, 0.905), fontsize=10, textcoords='axes fraction', annotation_clip=False, rotation=9)
ax[0].annotate('$\\tau^* \\ll 1$', xy=(0.25, 0.325), fontsize=10, textcoords='axes fraction', annotation_clip=False, rotation=45)

#ax[0].arrow(0.9, 0.8, -0.1, 0.2, overhang=0.4, head_width=0.025, head_length=0.05, lw=0.5, ec='k', fc='k', zorder=100, clip_on=False)

ax[0].arrow(0.25, 0.5, -0.14, 0.18, overhang=0.4, head_width=0.025, head_length=0.05, lw=0.5, ec='k', fc='k', zorder=100, clip_on=False)
#ax[0].annotate('$\\tau^*$', xy=(0.11, 0.72), fontsize=10, textcoords='axes fraction', annotation_clip=False, rotation=0)
ax[0].annotate('$\\tau^*$', xy=(0.25, 0.446), fontsize=12, textcoords='axes fraction', annotation_clip=False, rotation=0)
#ax[0].annotate('$\\tau^*$', xy=(0.175, 0.61), fontsize=10, textcoords='axes fraction', annotation_clip=False, rotation=0)





ax[1].arrow(0.18, 0.5, -0.12*0.6, 0.35*0.6, overhang=0.4, head_width=0.025, head_length=0.05, lw=0.5, ec='k', fc='k', zorder=100, clip_on=False)

#ax[1].arrow(0.5, 0.72, -0.1, 0.18, overhang=0.4, head_width=0.025, head_length=0.05, lw=0.75, ec='k', fc='k', zorder=100, clip_on=False)

ax[1].annotate('$\\nu \ll 1$', xy=(0.25, 0.30), fontsize=10, textcoords='axes fraction', annotation_clip=False, rotation=45)
ax[1].annotate('$\\nu$', xy=(0.17, 0.441), fontsize=12, textcoords='axes fraction', annotation_clip=False, rotation=0)



ax2 = ax[1].twinx()

ticks = [0,1]
ticklabels = ['$\displaystyle{\\frac{p_f^*}{\sigma_T^*} - 1}$', '$\displaystyle{\\frac{p_f^*}{\sigma_T^*}}$']

ax2.set_ylim([0,1])
ax2.set_yticks(ticks)
ax2.set_yticklabels(ticklabels)
ax2.set_frame_on(False)
plt.tight_layout()
#ax[1].set_aspect(1.0/ax[1].get_data_ratio(), adjustable='box')
plt.show()

'''
K=1.
period = sigma_T/(G+K)
time = np.linspace(0, period, 10001)
plt.plot(time/period, periodic_solution(time, G, K))
#plt.plot(time, -1+time*(G+K)/(1.+K))
plt.ylim([-1,0])
'''
