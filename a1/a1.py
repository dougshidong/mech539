#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def initCond(x):
    return (1.0 + np.tanh(250.0 * (x - 20.0))) / 2.0
def exactSol(x, c, t):
    return initCond(x - c * t)

# Space Discretization
L = 40.0
nx = 41
dx = L / (nx - 1)
x = np.linspace(0.0, L + dx, nx, dtype = float)

# Time Discretization
tf = 10.0
dt = 1.0
nt = int(np.ceil(tf / dt))

# Initial Condition
u0 = initCond(x)
# Constant
c = 0.5;

def upwind(u, dt, dx):
    u1 = list(u)
    for i in range(1, len(u) - 1):
        u1[i] = u[i] - dt / dx * c * (u[i] - u[i - 1])
    return u1

def lax(u, dt, dx):
    u1 = list(u)
    for i in range(1, len(u) - 1):
        u1[i] = ( u[i + 1] + u[i - 1] \
                - c * dt / dx * (u[i + 1] - u[i - 1]) ) / 2
    return u1

def laxwendroff(u, dt, dx):
    u1 = list(u)
    for i in range(1, len(u) - 1):
        u1[i] = u[i] - c * dt / dx * (u[i + 1] - u[i - 1]) / 2 \
                     + (c * dt / dx)**2 * (u[i + 1] - 2 * u[i] + u[i - 1]) / 2
    return u1

def leapfrog(u, u_old, dt, dx):
    u1 = list(u)
    for i in range(1, len(u) - 1):
        u1[i] = u_old[i] - c * dt / dx * (u[i + 1] - u[i - 1])
               
    return u1

def maccormack(u, dt, dx):
    u1 = list(u)
    u2 = list(u)
    for i in range(1, len(u) - 1):
        u1[i] = u[i] - c * dt / dx * (u[i + 1] - u[i])
    for i in range(1, len(u) - 1):
        u2[i] = ( u[i] + u1[i] - c * dt / dx * (u1[i] - u1[i - 1]) ) / 2

    return u2

# Scheme
# 0 = Upwind
# 1 = Lax
# 2 = Lax-Wendroff
# 3 = Leap-Frog
# 4 = MacCormack
# map the inputs to the function blocks
nscheme = 6
scheme = {1 : upwind,
          2 : lax,
          3 : laxwendroff,
          4 : leapfrog,
          5 : maccormack}

usol = np.empty([nscheme, nx])
# Exact Solution
usol[0, :] = exactSol(x, c, tf)
# Numerical Solutions
for s in range(1, nscheme):
    if s != 4:
        u = list(u0)
        for t in range(nt):
            print t
            u = scheme[s](u, dt, dx)
    else:
        u = list(u0)
        u_old = list(u0)
        for t in range(nt):
            u_new = scheme[s](u, u_old, dt, dx)
            u_old = u
            u = u_new
    usol[s, :] = u

# Plotting Q1
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(figsize=(10,5))
col=('b', 'g', 'r', 'c', 'm', 'y', 'k')
sname = ('Exact', 'Upwind', 'Lax', 'Lax-Wendroff', 'Leap-Frog', 'MacCormack')
for s in range(nscheme):
    plt.plot(x, usol[s, :], color = col[s],
             marker = 'o', mec = col[s], mfc = 'None', ms = 3,
             label = sname[s])

plt.axis([0, 40, -0.4, 1.1])

plt.title(r'Velocity Distribution, $\Delta x = 1.0$, $\Delta t = 1.0$')
plt.xlabel(r'$x$')
plt.ylabel(r'$u$')
plt.legend(loc='upper left')

plt.savefig('plot.pdf')
#plt.savefig('./report/Figures/q2.pdf')

# Plotting Q2
