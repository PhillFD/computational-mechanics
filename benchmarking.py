import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from model import *

def plot_benchmark():
    # Numerical ODE solution
    t_num, p_num = solve_ode(ode_model, 0.0, 10.0, 0.1, 0, [-1,1,1,0])

    # Analytical ODE solution
    t_analy = np.arange(0,10.1,0.1)
    p_analy = 1-np.e**(-1*t_analy)

    # Error Analysis
    err = []
    for i in range(101):
        err.append(abs(p_analy[i]-p_num[i])/p_analy[i])

    # Convergence Testing
    p_at_ten = []
    inverse_step_size = np.arange(1.0,3.1,0.1)
    step_size = 1/inverse_step_size

    for step in step_size:
        t_range, numerical_x = solve_ode(ode_model, 0, 10.0, step, 0, [-1,1,1,0])
        p_at_ten.append(numerical_x[-1])

    # Benchmark plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(12,6))
    ax1.plot(t_num, p_num, "-x", color="blue", label="numerical solution")
    ax1.plot(t_analy, p_analy, "-", color="red", label="analytical solution")
    ax1.set_title("benchmark: a=-1.00,b=1.00,P0=0,q0=1.00")
    ax1.set_xlabel("$t$")
    ax1.set_ylabel("$p$")
    ax1.legend()

    # Error analysis plot
    ax2.plot(t_num, err, "o", color="black")
    ax2.set_title("error analysis")
    ax2.set_xlabel("$t$")
    ax2.set_ylabel("relative error against benchmark")

    # Convergence test plot
    ax3.plot(inverse_step_size, p_at_ten, "x", color="black")
    ax3.set_title("timestep convergence")
    ax3.set_xlabel("1/∆$t$")
    ax3.set_ylabel("$P$($t$ = 10)")
    plt.show()