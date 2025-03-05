import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.loadtxt('results/N0001/average1D_VpRhoTS.txt')

# Extract columns 2 and 3
y = data[:, 1]  # Column 2
x = data[:, 2]  # Column 3

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(x, y, color='tab:blue')
plt.xlabel('Sound speed (m/s)')
plt.ylabel('Depth (m)')
plt.title('Sound speed profile at H03')
plt.ylim([5000,0])
plt.savefig('H03_speed.png',dpi=400)

