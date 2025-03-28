import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_convergence(adaptive_conv_filename, unadaptive_conv_filename):
    adap_data = pd.read_csv(adaptive_conv_filename)
    
    unadap_data = pd.read_csv(unadaptive_conv_filename)

    plt.plot(adap_data['Iteration'], adap_data['Proximity'], label='Adaptive Parallel Tempering')
    plt.plot(unadap_data['Iteration'], unadap_data['Proximity'], label='Vanilla Parallel Tempering')
    
    plt.xlabel("Iteration")
    plt.ylabel("Proximity to Optimum")
    plt.title("Solution Proximity")
    plt.legend()
    plt.grid(True)    
    plt.savefig('Solution_proximity.png')
    plt.show()

plot_convergence("adap_proximity.csv", "unadap_proximity.csv")