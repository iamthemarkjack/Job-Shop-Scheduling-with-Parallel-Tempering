import numpy as np
import matplotlib.pyplot as plt

# Reading the history
arr = []
fname = r'history.txt'
with open(fname) as f:
    f.readline() # Reading the header
    for line in f:
        arr.append(list(line.split(',')))
arr = np.array(arr, dtype = float)

no_replicas = max(arr[:,1].astype(int)) + 1 # No. of replicas
temp_his = [[] for _ in range(no_replicas)]

# Storing the temperature history in temp_his
for entry in arr:
    replica_idx = entry[1].astype(int)
    temp_his[replica_idx].append(entry[2])

# Plotting the temperature history
plt.figure(figsize = (12,8))
for i, his in enumerate(temp_his):
    plt.plot(his, label = f'Replica - {i}')
plt.legend()
plt.xlabel('Iterations')
plt.ylabel('Temperature')
plt.title('Temperature vs Iterations for each replica')
plt.ylim([0.5,7])
plt.savefig('Temp_hist.png')
plt.show()