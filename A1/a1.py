#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def initCond(x):
    return (1.0 + np.tanh(250.0 - (x - 20.0))) / 2.0

L = 40.0
nx = 41
dx = (nx - 1) / L

tf = 10.0
dt = 1.0

x = np.arange(0.0, L, dx, dtype = float)
u0 = initCond(x);

plt.figure()
plt.plot(u0)
