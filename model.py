import numpy as np
from matplotlib import pyplot as plt

def ode_model(t, p, q, a, b, p0):
    # Reservoir pressure ODE
    return -a*q-b*(p-p0)

def solve_ode(f, t0, t1, dt, p0, pars):
    # Improved euler method solver
    t = np.arange(t0,t1+dt,dt)
    p = 0.*t
    p[-1] = p0
    for i in range(0, len(t)):
        f0=f(t[i-1],p[i-1],*pars)
        p[i] = p[i-1] + dt*(f0+f(t[i],p[i-1] + dt*f0, *pars))/2
    return t, p

def load_reservoir_pressure():
    # Load time and pressure from file
    t, Pa = np.genfromtxt('gs_pres.txt', delimiter=',', skip_header=1).T
    return t, Pa

def load_reservoir_massflow():
     # Load time and mass flow from file
    t, q = np.genfromtxt('gs_mass.txt', delimiter=',', skip_header=1).T
    return t, q
