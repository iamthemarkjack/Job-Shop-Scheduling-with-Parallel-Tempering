import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_parallel_tempering_history(sap_filename, temp_filename):
    # Load the SAP history file
    sap_data = pd.read_csv(sap_filename)
    temp_data = pd.read_csv(temp_filename)
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot SAP History
    iterations = sap_data.iloc[:, 0]  # First column is iteration numbers
    for col in sap_data.columns[1:]:  # Remaining columns are SAP values for chain pairs
        ax1.plot(iterations, sap_data[col], label=col)
    
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Swap Acceptance Probability")
    ax1.set_title("Swap Acceptance Probability History")
    ax1.legend()
    ax1.grid(True)
    
    # Plot Temperature History
    iterations = temp_data.iloc[:, 0]  # First column is iteration numbers
    for col in temp_data.columns[1:]:  # Remaining columns are temperature values for chains
        ax2.plot(iterations, temp_data[col], label=col)
    
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Temperature")
    ax2.set_title("Temperature History")
    ax2.legend()
    ax2.grid(True)
    
    # Adjust layout and display
    plt.tight_layout()
    plt.show()

plot_parallel_tempering_history("sap_history.csv", "temp_history.csv")