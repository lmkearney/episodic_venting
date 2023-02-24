import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.integrate import trapz

rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
rc("text", usetex=True)

sand = (0.9, 0.85, 0.0)
mud = (0.66, 0.41, 0.21)
p_f = 0.5
sigma_T = 0.25


def MandelRoots(N):
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

def p_m_transition(x, t, N):
    p = 0.5 * x * (2.0 - x)
    for i in range(1, N):
        p += (
            16
            / np.pi**3
            * (-1) ** i
            / (2 * i - 1.0) ** 3
            * np.cos((2 * i - 1) * np.pi / 2.0 * (1.0 - x))
            * (np.exp(-0.25 * np.pi**2 * (2 * i - 1.0) ** 2 * t))
        )
    return p


def p_m_avg_transition(t, N):
    p = 1 / 3.0
    for i in range(1, N):
        p -= (
            32
            / np.pi**4
            * np.exp(-t * (2.0 * i - 1) ** 2 * np.pi**2 / 4.0)
            / (2.0 * i - 1.0) ** 4
        )
    return p

J = 101
x_grid = np.linspace(0, 1, J)
dx = x_grid[1] - x_grid[0]

lnw = 1.25


N = 100
G = 2.0
K = 5.0

Ts = 10000
t_grid = np.linspace(0, 2.8, Ts)
dt = t_grid[1] - t_grid[0]

sigma_T = 0.25
p_f = 0.5

alpha = dt / (2.0 * dx**2)
beta = K * dt / (2.0 * dx)

# initial condition
U = np.zeros(J)

# matrix elements
A_u = (
    np.diagflat([-alpha for i in range(J - 1)], -1)
    + np.diagflat(
        [(1 + 2 * alpha) + alpha / beta * (1 - beta)]
        + [1.0 + 2.0 * alpha for i in range(J - 2)]
        + [1.0 + alpha]
    )
    + np.diagflat([-alpha] + [-alpha for i in range(J - 2)], 1)
)

B_u = (
    np.diagflat([alpha for i in range(J - 1)], -1)
    + np.diagflat(
        [(1 - 2 * alpha) + alpha / beta * (1 + beta)]
        + [1.0 - 2.0 * alpha for i in range(J - 2)]
        + [1.0 - alpha]
    )
    + np.diagflat([alpha] + [alpha for i in range(J - 2)], 1)
)

C_t = np.zeros(J) + 1.0
C_t[0] += G * alpha / beta
U_record = []
U_record.append(U)

avg_pm = np.zeros(len(t_grid))
ps = np.zeros(len(t_grid))

t_stop = 1.4

for ti in range(1, Ts):
    U_new = np.linalg.solve(A_u, B_u.dot(U) + dt * C_t)
    if U_new[0] >= p_f:
        U_new[0] = p_f - sigma_T
    U = U_new
    U_record.append(U)


for i in range(len(U_record)):
    avg_pm[i] = trapz(U_record[i], x_grid)
    ps[i] = U_record[i][0]

P_s = ps
P_m = avg_pm
fig, ax = plt.subplots(1, 1, figsize=(8, 3))

ax.axhline(p_f, lw=0.75, c="k", alpha=1)
ax.axhline(p_f - sigma_T, lw=0.75, c="k", alpha=1)

ax.plot(t_grid, P_s, c="k", alpha=0.15, lw=1.0, label="sandstone")
ax.plot(t_grid, P_m, c="k", alpha=0.4, lw=1.0, label="mudstone", ls=(0, (5, 0.75)))

t1 = t_grid + 0.3917
ax.plot(t1, 0.4147 + p_m_avg_transition(t_grid, N), c=mud, ls=(0, (5, 0.75)), lw=1.25)
ax.plot(t1, np.full_like(t1, 0.4147), c=sand, lw=1.25)
ax.set_xlabel("time $t^*$")
ax.set_ylabel(r"average overpressure $\overline{p^*}$")

ax.set_xlim([0.0, t_grid[-1]])
ax.set_ylim([0.0, 1.0])

ax2 = ax.twinx()

ticks = [0.25, 0.5]
ticklabels = [r"$p_f^* - \sigma_T^*$", "$p_f^*$"]

ax2.set_yticks(ticks)
ax2.set_yticklabels(ticklabels)
ax2.set_frame_on(False)
plt.tight_layout()
# plt.savefig('./sample.pdf', bbox_inches="tight")
plt.show()
