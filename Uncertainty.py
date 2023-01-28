import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import multivariate_normal
from Calibration import *

def Uncertainty():
    #get best fit
    init_pars = [-4.59528059079902e-06, 3.480663889651829]  # optimum parameters
    p0 = 25.16
    tP, pa1 = load_reservoir_pressure()
    tM, pa2 = solve_reservoir_ode(ode_model, 2009, 2019, 0.01, p0, init_pars)

    f, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.plot(tP, pa1, 'b-', label='observations')  # observations
    ax.plot(tM, pa2, 'g-', label='best fit')  # initial model
    ax.set_xlabel('time, $t$ [years]')
    ax.set_ylabel('pressure, $P$ [MPa]]')
    ax.legend()
    plt.show()

    sigma = [0.25]*len(tM)
    x = np.linspace(np.min(tM), np.max(tM), 101)
    p, cov = curve_fit(solve_ode, tM, pa2, sigma=sigma)
    ps = np.random.multivariate_normal(p, cov, 100)  # mutliple variances

    for pi in ps:
        ax.plot(x, ode_model(x, *pi), 'k-', alpha=0.2, lw=0.5)
    ax.plot([], [], 'k-', lw=0.5, label='posterior samples')
    ax.legend()

if __name__ == '__main__':
    Uncertainty()