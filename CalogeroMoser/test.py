import numpy as np
import matplotlib.pyplot as plt

zeros = np.load("Output/-1.00.npy")

X = zeros[:, 0]
Y = zeros[:, 1]

plt.scatter(X,Y)
plt.show()
