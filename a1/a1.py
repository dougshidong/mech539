#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def initCond(x):
    return (1.0 + np.tanh(250.0 * (x - 20.0))) / 2.0

def exactSol(x, c, t):
    return initCond(x - c * t)

def solError(un, ue):
    return np.linalg.norm(un - ue) / len(un)**0.5

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

# Space Discretization
L = 40.0
nx = 41
dx = L / nx
# Time Discretization
tf = 10.0
dt = 1.0
nt = int(np.ceil(tf / dt))
# Constant
c = 0.5;

q = 4
q3calc = 0
# Exact Solution
plt.figure(figsize=(10,5))
if q == 1 or q == 2:
    xe = np.linspace(0.0, L, 99500, dtype = float)
    uexact = exactSol(xe, c, tf)
    plt.plot(xe, uexact, color = 'k', label = 'Exact')

CFL = 0.5

if q == 1:
    nscheme = 5
    nxlist = [41]
    scheme = {0 : upwind,
              1 : lax,
              2 : laxwendroff,
              3 : leapfrog,
              4 : maccormack}
if q == 2:
    nscheme = 2
    nxlist = [41, 81, 161, 321]
    scheme = {0 : upwind,
              1 : maccormack}
if q == 3:
    nscheme = 5
    nxlist = [41, 81, 161, 321, 641, 1281, 2561, 5121, 10241, 20481, 40561]
    err = np.empty((nscheme, len(nxlist)))
    scheme = {0 : upwind,
              1 : lax,
              2 : laxwendroff,
              3 : leapfrog,
              4 : maccormack}
if q == 4:
    nscheme = 1
    nxlist = [41, 81, 161, 162, 167]
    #scheme = {0 : laxwendroff}
    scheme = {0 : lax}

for (nxi, nx) in enumerate(nxlist):
    dx = L / (nx - 1)
    dt = CFL * dx / c
    if(q == 4):
        dt = 0.5
        CFL = c * dt / dx
        print 'CFL: ', CFL
    nt = int(np.ceil(tf / dt))

    usol = np.empty((nscheme, nx))
    x = np.linspace(0.0, L, nx, dtype = float)
    # Initial Condition
    u0 = initCond(x)
    
    if(q != 3 or q3calc == 1):
        # Numerical Solutions
        for s in range(0, nscheme):
            if s != 3:
                u = list(u0)
                for t in range(nt):
                    u = scheme[s](u, dt, dx)
            else:
                u_old = list(u0)
                u = scheme[0](u_old, dt, dx)
                for t in range(1, nt):
                    u_new = scheme[s](u, u_old, dt, dx)
                    u_old = u
                    u = u_new
            usol[s, :] = u
    
    # Plotting
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    col=('b', 'g', 'r', 'c', 'm', 'y', 'k')
    symb = ('o', 'v', '<', 's', 'p', '*', 'h', 'H', 'D', 'd')

    if q == 1:
        sname = ('Upwind', 'Lax', 'Lax-Wendroff', 'Leap-Frog', 'MacCormack')
        for s in range(nscheme):
            plt.plot(x, usol[s, :], color = col[s],
                     marker = symb[s], mec = col[s], mfc = 'None', ms = 6,
                     label = sname[s])
        plt.axis([0, 40, -0.4, 1.1])
        plt.title(r'Velocity Distribution, $\Delta x = %.1f$, $\Delta t = %.1f$' % (dx, dt))
        plt.legend(loc='upper left')
        plt.xlabel(r'$x$')
        plt.ylabel(r'$u$')

    if q == 2:
        sname = ('Upwind', 'MacCormack')
        for s in range(nscheme):
            plt.plot(x, usol[s, :], color = col[nxi],
                     marker = symb[s], mec = col[nxi], mfc = 'None', ms = 6,
                     label = sname[s] + ' nx = ' + str(nx))
        plt.axis([20, 30, -0.4, 1.1])
        plt.title(r'Grid Study Velocity Distribution')
        plt.legend(loc='upper left', fontsize = 10)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$u$')
    if q == 3:
        if q3calc == 1:
            print nxlist[nxi]
            uexact = exactSol(x, c, tf)
            for s in range(0, nscheme):
                err[s, nxi] = solError(usol[s, :], uexact)
    if q == 4:
        #sname = (['Lax-Wendroff',''])
        sname = (['Lax',''])
        for s in range(nscheme):
            plt.plot(x, usol[s, :], color = col[nxi],
                     marker = symb[s], mec = col[nxi], mfc = 'None', ms = 6,
                     label = sname[s] + ' CFL = ' + str(CFL))
        plt.axis([15, 30, -0.4, 1.1])
        plt.title(r'Stability Condition')
        plt.legend(loc='upper left', fontsize = 10)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$u$')

if q == 3:
    if q3calc == 1:
        np.save('q3resultscfl08', err)
    else:
        err = np.load('q3resultscfl08.npy')
    sname = ('Upwind', 'Lax', 'Lax-Wendroff', 'Leap-Frog', 'MacCormack')
    for s in range(nscheme):
        plt.loglog(nxlist, err[s, :], color = col[s],
                 marker = symb[s], mec = col[s], mfc = 'None', ms = 6,
                 label = sname[s])
    plt.axis([0, 42000, 0, 0.2])
    plt.title(r'Order of Accuracy')
    plt.legend(loc='upper right')
    plt.xlabel(r'Number of Points')
    plt.ylabel(r'Error')

plt.savefig('plot.pdf')
plt.savefig('./report/Figures/q4_1.pdf')

# Plotting Q2
