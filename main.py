from model import *
from Calibration import *
from prediction import *
from benchmarking import *

def main(): 
    # plot_scenario_ode() takes the rate of operation as input (i.e. 2 as double the operation rate)
    # default is 1 - no change in operation rate 

    plot_benchmark()
    plot_curve_fit()
    plot_reservoir_model()
    plot_bestfit()
    plot_scenario_ode(1) # predicts pessure over next 10 years with varying operation rates
    plot_all_scenarios()
    plot_all_leakage_vs_pressure()

if __name__ == "__main__":
    main()