import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_parallel_tempering_history(sap_filename, temp_filename):
    sap_data = pd.read_csv(sap_filename)
    temp_data = pd.read_csv(temp_filename)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    iterations = sap_data.iloc[:, 0]
    for col in sap_data.columns[1:]:
        ax1.plot(iterations, sap_data[col], label=col)
    
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Swap Acceptance Probability")
    ax1.set_title("Swap Acceptance Probability History")
    ax1.legend()
    ax1.grid(True)
    
    iterations = temp_data.iloc[:, 0]
    for col in temp_data.columns[1:]:
        ax2.plot(iterations, temp_data[col], label=col)
    
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Temperature")
    ax2.set_title("Temperature History")
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('SAP_Convergence.png')
    plt.show()

plot_parallel_tempering_history("adap_sap_history.csv", "adap_temp_history.csv")