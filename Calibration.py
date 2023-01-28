from scipy.optimize import curve_fit
from model import *

def interpolate_reservoir_massflow(t):
    # Interpolates q (massflow) against an arbitrary time interval
    tM, q = load_reservoir_massflow()
    return np.interp(t, tM, q)

def solve_reservoir_ode(f, t0, t1, dt, p0, pars):
    # Original Improved Euler Method taking into account q (massflow)
    t = np.arange(t0, t1+dt, dt)
    p = 0.*t
    p[-1] = p0
    q = interpolate_reservoir_massflow(t)
    for i in range(0,len(t)):
        f0 = f(t[i-1], p[i-1], q[i-1], *pars, p0)
        p1 = p[i-1] + dt*f0
        f1 = f(t[i], p1, q[i], *pars, p0)
        p[i] = p[i-1] + dt*(f0+f1)/2
    return t, p

def solve_curve_fit(f, t, p0, pars):
    # Improved Euler Method for Curve Fitting
    dt = t[1] - t[0]
    p = 0.*t
    p[-1] = p0
    q = interpolate_reservoir_massflow(t)
    for i in range(0,len(t)):
        f0 = f(t[i-1], p[i-1], q[i-1], *pars, p0)
        p1 = p[i-1] + dt*f0
        f1 = f(t[i], p1, q[i], *pars, p0)
        p[i] = p[i-1] + dt*(f0+f1)/2
    return t, p

def plot_curve_fit():
    # Returns optimum parameters for the best fitting model to our data for a given prediction
    to,Po = load_reservoir_pressure()

    # optimum parameters [-4.59528059079902e-06,3.480663889651829]
    # our guess
    pars=[-7.55e-6,4.95]
    p0=25.16
    tm,Pm = solve_curve_fit(ode_model, to, p0, pars)

    def Tmodel(t, *pars):
        p0=25.16
        tm,Pm = solve_curve_fit(ode_model, t, p0, pars)
        return Pm

    # our guess
    p1=[-7.55e-6,4.95]
    constants=curve_fit(Tmodel, to, Po, p1)
    a_const=constants[0][0]
    b_const=constants[0][1]
    print(a_const)
    print(b_const)
    pars=[a_const, b_const]
    p0=25.16
    tmi,Pmi = solve_curve_fit(ode_model, to, p0, pars)

    f,ax = plt.subplots(1, 1, figsize=(12,6))
    ax.plot(to,Po, 'ko', label='observations')
    ax.plot(tm,Pm, 'r-', label='model guess')
    ax.plot(tmi,Pmi, 'b-', label='model improved')
    ax.set_xlabel('time, $t$ [years]')
    ax.set_ylabel('pressure, $P$ [MPa]')
    ax.legend()
    plt.show()

def plot_reservoir_model():
    # Plotting the reservoir ODE with varying q 
    init_pars = [-7.55e-6,4.95] # ad-hoc parameters 
    p0=25.16
    tP, pa1 = load_reservoir_pressure()
    tM, pa2 = solve_reservoir_ode(ode_model,2009,2019,0.25,p0,init_pars) 
    pressure_misfit = pa2-pa1

    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,6))

    # plot pressure vs time
    ax1.plot(tP, pa1, 'b-', label='observations') # observations
    ax1.plot(tM, pa2, 'g-', label='a = -7.55e-6, b = 4.95') # initial model
    ax1.set_xlabel('time, $t$ [years]')
    ax1.set_ylabel('pressure, $P$ [MPa]]')
    ax1.set_title('best fit LPM model')
    ax1.legend()

    # plot pressure misfit 
    ax2.plot(tP, pressure_misfit, 'rX') # initial model
    ax2.axhline(y=0, linestyle='dashed', color='black')
    ax2.set_xlabel('time, $t$ [years]')
    ax2.set_ylabel('pressure misfit, $P$ [MPa]]')
    ax2.set_title('best fit LPM model')
    plt.show()

def plot_bestfit():
    # Plotting the reservoir ODE with varying q 
    init_pars = [-4.59528059079902e-06,3.480663889651829] # optimum parameters 
    p0=25.16
    tP, pa1 = load_reservoir_pressure()
    tM, pa2 = solve_reservoir_ode(ode_model,2009,2019,0.25,p0,init_pars) 
    pressure_misfit = pa2-pa1

    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,6))

    # plot pressure vs time
    ax1.plot(tP, pa1, 'b-', label='observations') # observations
    ax1.plot(tM, pa2, 'g-', label='a = -4.59e-06, b = 3.48') # best fit 
    ax1.set_xlabel('time, $t$ [years]')
    ax1.set_ylabel('pressure, $P$ [MPa]]')
    ax1.set_title('best fit LPM model')
    ax1.legend()

    # plot pressure misfit 
    ax2.plot(tP, pressure_misfit, 'rX') # best fit 
    ax2.axhline(y=0, linestyle='dashed', color='black')
    ax2.set_xlabel('time, $t$ [years]')
    ax2.set_ylabel('pressure misfit, $P$ [MPa]]')
    ax2.set_title('best fit LPM model')
    plt.show()