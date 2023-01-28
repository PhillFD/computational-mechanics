from Calibration import *

def solve_scenario_ode(f, t0, t1, dt, p0, pars, rate=1):
    # Check for valid rate 
    if rate < 0:
        raise ValueError("Rate must be positive")
    if rate > 100:
        raise ValueError("Rate must be realistic")
   
    # Original Improved Euler Method taking into account different operation rates of q (massflow)
    t = np.arange(t0, t1+dt, dt)
    qt = np.arange(t0, 2019+dt, dt)
    p = 0.*t
    p[-1] = p0
    q1 = interpolate_reservoir_massflow(qt) # interpolating mass flow for first 10 years
    q2 = interpolate_reservoir_massflow(qt) * rate # interpolating mass flow for next 10 years
    q = np.concatenate([q1,q2])
    for i in range(0,len(t)):
        f0 = f(t[i-1], p[i-1], q[i-1], *pars, p0)
        p1 = p[i-1] + dt*f0
        f1 = f(t[i], p1, q[i], *pars, p0)
        p[i] = p[i-1] + dt*(f0+f1)/2
    return t, p

def plot_scenario_ode(rate=1):
    # Plotting the reservoir model prediction for the next 10 years at varying operation capacities
    opt_pars=[-4.59528059079902e-06,3.480663889651829] # optimum parameters
    p0=25.16
    t1, pa1 = load_reservoir_pressure()
    t2, pa2 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,rate) # best fit model for next 10 years

    f,ax = plt.subplots(1, 1, figsize=(12,6))
    ax.plot(t1, pa1, 'b-', label='observations') # observations
    ax.plot(t2, pa2, 'r-', label='model prediction') # best fit model 
    ax.set_xlabel('time, $t$ [years]')
    ax.set_ylabel('pressure, $P$ [MPa]]')
    ax.set_title('Model prediction for the next 10 years at ' + str(rate) + 'x operation capacity')
    ax.legend()
    plt.show()

def plot_all_scenarios():
    # Plotting all reservoir model predictions for the next 10 years at varying operation capacities
    opt_pars=[-4.59528059079902e-06,3.480663889651829] # optimum parameters
    p0=25.16
    t1, pa1 = load_reservoir_pressure()
    t2, pa2 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars) # no change in operation capacity
    t3, pa3 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,2) # double in operation capacity
    t4, pa4 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,0.05) # decrease in operation capacity
    t5, pa5 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,1.25) # optimal operation capacity

    t3 = np.ma.masked_less(t3, 2019) # only show predicted data
    t4 = np.ma.masked_less(t4, 2019) # hide data before 2019
    t5 = np.ma.masked_less(t5, 2019) 

    f,ax = plt.subplots(1, 1, figsize=(12,6))
    ax.plot(t1, pa1, 'b-', label='observations') 
    ax.plot(t2, pa2, 'r-', label='standard operation capacity') 
    ax.plot(t3, pa3, 'g-', label='double operation capacity') 
    ax.plot(t4, pa4, 'y-', label='decreased operation capacity') 
    ax.plot(t5, pa5, 'm-', label='optimal operation capacity') 
    plt.xlabel('time, $t$ [years]')
    plt.ylabel('pressure, $P$ [MPa]]')
    plt.title('Ahuroa LPM: what-if scenarios')
    plt.legend()
    plt.show()

def plot_all_leakage_vs_pressure():
    # Plotting all reservoir model predictions for the next 10 years at varying operation capacities
    opt_pars=[-4.59528059079902e-06,3.480663889651829] # optimum parameters
    p0=25.16

    t1, pa1 = load_reservoir_pressure()
    t2, pa2 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars) # no change in operation capacity
    t3, pa3 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,2) # double in operation capacity
    t4, pa4 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,0.05) # decrease in operation capacity
    t5, pa5 = solve_scenario_ode(ode_model,2009,2029,0.25,p0,opt_pars,1.25) # optimal operation capacity

    t2 = np.ma.masked_less(t2, 2019) # only show predicted data
    t3 = np.ma.masked_less(t3, 2019) # hide data before 2019
    t4 = np.ma.masked_less(t4, 2019) 
    t5 = np.ma.masked_less(t5, 2019) 

    # Plot Pressure Leakage against time 
    threshold_pressure = 25.16169799530975 # threshold pressure = p0
    k = 3.480663889651829/threshold_pressure # leakage co-efficient
    threshold_over_pressure = 0. # overpressure value (p-p0)

    overpressure_array_2 = pa2 - threshold_pressure # array of overpressure values normal
    pressure_leakage_2 = 0.*overpressure_array_2

    overpressure_array_3 = pa3 - threshold_pressure # array of overpressure values double
    pressure_leakage_3 = 0.*overpressure_array_3

    overpressure_array_4 = pa4 - threshold_pressure # array of overpressure values decrease
    pressure_leakage_4 = 0.*overpressure_array_4

    overpressure_array_5 = pa5 - threshold_pressure # array of overpressure values optimal
    pressure_leakage_5 = 0.*overpressure_array_5

    # if pressure leakage goes over threshold overpressure, then calculate pressure leakage
    for i in range(len(overpressure_array_2)):
        if overpressure_array_2[i] > threshold_over_pressure:
            pressure_leakage_2[i] = k*(overpressure_array_2[i])
        else:
            pressure_leakage_2[i] = 0    

    for i in range(len(overpressure_array_3)):
        if overpressure_array_3[i] > threshold_over_pressure:
            pressure_leakage_3[i] = k*(overpressure_array_3[i])     
        else:
            pressure_leakage_3[i] = 0

    for i in range(len(overpressure_array_4)):
        if overpressure_array_4[i] > threshold_over_pressure:
            pressure_leakage_4[i] = k*(overpressure_array_4[i])
        else:
            pressure_leakage_4[i] = 0

    for i in range(len(overpressure_array_5)):
        if overpressure_array_5[i] > threshold_over_pressure:
            pressure_leakage_5[i] = k*(overpressure_array_5[i])  
        else:
            pressure_leakage_5[i] = 0
        
    fig, ax1 = plt.subplots(1,1,figsize=(12,6))
    ax1.plot(t2, pressure_leakage_2, 'b-', label='normal capacity')
    ax1.plot(t3, pressure_leakage_3, 'r-', label='double capacity')
    ax1.plot(t4, pressure_leakage_4, 'g-', label='decreased capacity')
    ax1.plot(t4, pressure_leakage_5, 'm-', label='optimal capacity')
    ax1.set_xlabel('time, $t$ [years]')
    ax1.set_ylabel('pressure leakage, $P$ [MPa]]')
    ax1.set_title('Ahuroa LPM: what-if scenarios')
    ax1.legend()
    plt.show()

  

