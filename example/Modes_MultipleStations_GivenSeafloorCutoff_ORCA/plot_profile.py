import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.loadtxt('results/H03N1/average1D_VpRhoTS.txt')

# Extract columns 1 and 2
y = data[:, 0]  # Column 1 (Depth)
x = data[:, 1]  # Column 2 (Sound speed)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(x, y, color='tab:blue')
plt.xlabel('Sound speed (m/s)')
plt.ylabel('Depth (m)')
plt.title('Sound speed profile at H03')

# Invert y-axis so depth increases downward
plt.gca().invert_yaxis()

plt.show()
# plt.savefig('H03_speed.png', dpi=400)

