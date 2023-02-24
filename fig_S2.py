import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as np
from scipy.integrate import trapz

rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
rc("text", usetex=True)

sand = (0.9, 0.85, 0.0)
mud = (0.66, 0.41, 0.21)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


J = 101
x_grid = np.linspace(0, 1, J)
dx = x_grid[1] - x_grid[0]

G = 2.0
K = 5.0
p_f = 0.5
sigma_T = 0.25

Ts = 10000
t_grid = np.linspace(0, 2.8, Ts)
dt = t_grid[1] - t_grid[0]

fig, ax = plt.subplots(1, 1, figsize=(8, 3))

# Crank-Nicolson simulation method
alpha = dt / (2.0 * dx**2)
beta = K * dt / (2.0 * dx)

# initial condition
U = np.zeros(J)  # +p_f

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
        # print ti*dt
    if ti * dt >= t_stop:
        C_t = np.zeros(J)
    U = U_new
    U_record.append(U)


for i, value in U_record:
    avg_pm[i] = trapz(value, x_grid)
    ps[i] = value[0]

P_s = ps
P_m = avg_pm

ax.axhline(p_f, lw=0.75, c="k", alpha=1)
ax.axhline(p_f - sigma_T, lw=0.75, c="k", alpha=1)
ax.axvline(t_stop, lw=1.0, c="k", ls="dotted", alpha=1)

ax.plot(t_grid, G * t_grid, alpha=0.15, c="k", lw=1.0)
ax.plot(t_grid, t_grid, alpha=0.3, c="k", lw=1.0, ls=(0, (5, 0.75)))
ax.plot(t_grid, P_s, c=sand, lw=1.0, label="sandstone")
ax.plot(t_grid, P_m, c=mud, lw=1.0, label="mudstone", ls=(0, (5, 0.75)))


# now compute pressure if compression continued

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

for ti in range(1, Ts):
    U_new = np.linalg.solve(A_u, B_u.dot(U) + dt * C_t)
    if U_new[0] >= p_f:
        U_new[0] = p_f - sigma_T
        # print ti*dt
    U = U_new
    U_record.append(U)


for i, value in enumerate(U_record):
    avg_pm[i] = trapz(value, x_grid)
    ps[i] = value[0]

P_s = ps
P_m = avg_pm

idx = find_nearest(t_grid, t_stop)
ax.plot(t_grid[idx:], P_m[idx:], alpha=0.3, c="k", lw=1.0, ls=(0, (5, 0.75)))


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
plt.show()
