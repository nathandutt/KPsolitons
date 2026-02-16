import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.loadtxt('evolution.csv', delimiter=',', skiprows=1)
t, x, y = data[:, 0], data[:, 1], data[:, 2]

# Plot the trajectory as a line
plt.figure(figsize=(10, 6))
plt.plot(x, y, '-o', label='Trajectory')
plt.title('Trajectory over Time')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show()
