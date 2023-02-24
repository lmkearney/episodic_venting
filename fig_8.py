import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np

rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
rc("text", usetex=True)


def roots(K, N):
    eps = 0.00001
    zeros = np.zeros(N)
    for i in range(1, N):
        a1 = (i - 1) * np.pi + np.pi / 2 + eps
        a2 = a1 + np.pi - eps
        for j in range(40):
            y1 = np.tan(a1) + a1 / K
            y2 = np.tan(a2) + a2 / K
            am = (a1 + a2) / 2.0
            ym = np.tan(am) + am / K
            if ym * y1 > 0:
                a1 = am
            else:
                a2 = am
            if abs(y2) < eps:
                am = a2
        zeros[i] = am
    return zeros


def p_s_sawtooth(t, G):
    ps = p_0 + np.zeros(len(t))
    dt = t[1] - t[0]
    for i in range(1, len(t)):
        ps[i] = ps[i - 1] + G * dt
        if ps[i] >= p_f:
            ps[i] = p_f - sigma_T
    return ps


def p_m_avg_comp(t, G, K):
    zeros = roots(K, N)
    p = t + (G - 1.0) * (
        (2 * t - 2 / 3.0) / (2.0 * (K + 1.0)) + K / (3.0 * (1.0 + K) ** 2)
    )
    for i in range(1, N):
        p -= (
            2.0
            * K
            * (G - 1.0)
            * np.tan(zeros[i])
            * (np.exp(-zeros[i] ** 2 * t))
            / (zeros[i] ** 3 * (K**2 + K + zeros[i] ** 2))
        )
    return p


def p_s_comp(t, G, K):
    zeros = roots(K, N)
    p = t + (G - 1.0) * (t / (K + 1.0) + K / (3.0 * (1.0 + K) ** 2))
    for i in range(1, N):
        p -= (
            2.0
            * K
            * (G - 1.0)
            * np.exp(-zeros[i] ** 2 * t)
            / (zeros[i] ** 2 * (K**2 + K + zeros[i] ** 2))
        )
    return p


def p_s_pipe(t, G):
    ps = np.zeros(len(t))
    dt = t[1] - t[0]
    for i in range(1, len(t)):
        ps[i] = ps[i - 1] + G * dt
        if ps[i] >= p_f:
            ps[i] = p_f - sigma_T
    return ps


def p_m_avg_pipe(t, t0, K):
    zeros = roots(K, N)
    p = -sigma_T / (1.0 + K) * np.heaviside(t - t0, 1)
    for i in range(1, N):
        p -= (
            2.0
            * K
            * sigma_T
            * np.tan(zeros[i])
            * np.exp(-zeros[i] ** 2 * abs(t - t0))
            / (zeros[i] * (K**2 + K + zeros[i] ** 2))
            * np.heaviside(t - t0, 1)
        )
    return p


def p_s_pipe(t, t0, K):
    zeros = roots(K, N)
    p = -sigma_T / (1.0 + K) * np.heaviside(t - t0, 1)
    for i in range(1, N):
        p -= (
            2.0
            * K
            * sigma_T
            * np.exp(-zeros[i] ** 2 * abs(t - t0))
            / (K**2 + K + zeros[i] ** 2)
            * np.heaviside(t - t0, 1)
        )
    return p


sand = (0.9, 0.85, 0.0)
mud = (0.66, 0.41, 0.21)

G_nu = 1.0 + 10 ** np.linspace(-3, 3, 10000)
Contr = 1.0 / (G_nu)

fig, ax = plt.subplots(figsize=(5.55, 3))

G = 2.0
K = 8.0
time = np.linspace(0, 2.0, 10001)
N = 1000

sigma_T = 1.0
p_0 = 0.9
p_f = 1.0

P_s = p_0 + p_s_comp(time, G, K)
P_m = p_0 + p_m_avg_comp(time, G, K)

t_event = 0
while any(P_s > p_f):
    t_event = time[np.argmax(P_s > p_f)]
    P_s += p_s_pipe(time, t_event, K)
    P_m += p_m_avg_pipe(time, t_event, K)

ax.axvline(2.0, c="k", alpha=0.5, lw=0.75, ls="dotted")
ax.semilogx(G_nu, Contr, c="orangered", label="compression")

ax.plot(1.0 + K / G, 1.0 / (1.0 + K / G), "o", zorder=10, clip_on=False, c="royalblue")
ax.plot(1, 1, "o", zorder=10, clip_on=False, c="seagreen")

ax.set_xlim([1e0, 1e2])
ax.set_ylim([0, 1])
ax.set_ylabel(r"period $\displaystyle \\frac{\\tau \Gamma_s}{\sigma_T}$")
ax.set_xlabel(
    r"venting frequency multiplier $\displaystyle 1 + \\frac{\\nu}{\gamma} = 1 + \\frac{h_m}{h_s} \, \\frac{\\alpha_m}{\\alpha_s} \, \\frac{1 + \lambda_s/\mu_s}{1 + \lambda_m/\mu_m}$"
)

ax2 = ax.twinx()
ticks = np.linspace(0, 1, 6)
ax2.set_yticks(ticks)
ax2.set_ylabel("pressure recharge contribution \n from compression")
ax2.set_frame_on(False)

# inset figures
axins = inset_axes(
    ax,
    width=1,
    height=0.7,
    loc="lower left",
    bbox_to_anchor=(0.2, 0.55, 0, 0),
    bbox_transform=ax.transAxes,
)
axins2 = inset_axes(
    ax,
    width=1,
    height=0.7,
    loc="lower left",
    bbox_to_anchor=(0.6, 0.15, 0, 0),
    bbox_transform=ax.transAxes,
)

axins.plot(time, p_s_sawtooth(time, G), c=sand)
axins2.plot(time, P_s, c=sand)

axins.plot(1, 0, "o", zorder=10, clip_on=False, c="seagreen", ms=4.0)
axins2.plot(1, 0, "o", zorder=10, clip_on=False, c="royalblue", ms=4.0)

axins.set_xticks([])
axins.set_yticks([])
axins.set_xlabel("time")
axins.set_ylabel("$p_s$")
axins.set_xlim([1, 2])
axins.set_ylim([0, 1])

axins2.set_xticks([])
axins2.set_yticks([])
axins2.set_xlabel("time")
axins2.set_ylabel("$p_s$")
axins2.set_xlim([1, 2])
axins2.set_ylim([0, 1])

plt.tight_layout()
plt.show()
