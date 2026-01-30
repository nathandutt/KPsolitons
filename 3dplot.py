import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.colors as colors

# load data
data = np.loadtxt("pole_evolution.csv", delimiter=",")

t_arr = data[:, 0]
zdata = data[:, 1:]

n_poles = zdata.shape[1] // 2
n_steps = len(t_arr)

# Combine real and imaginary parts
all_Z = np.zeros((n_steps, n_poles), dtype=complex)
for ti in range(n_steps):
    for zi in range(0, zdata.shape[1], 2):
        all_Z[ti, zi // 2] = zdata[ti, zi] + 1j * zdata[ti, zi + 1]

# Spatial grid
x_range = 200
x_s = np.linspace(-x_range, x_range, 800)

def field(x, Z_s):
    return np.real(np.sum(1 / (x - Z_s)))

# Precompute phi
phi = np.zeros((n_steps, x_s.size))
for i in range(n_steps):
    for j in range(x_s.size):
        phi[i, j] = field(x_s[j], all_Z[i])
y_s = np.linspace(-300, 300, 1200)
X, Y = np.meshgrid(x_s, y_s)
Z = phi
Z_x, Z_y = np.gradient(phi)
#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
fig, ax = plt.subplots()
#surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                      linewidth=0, antialiased=False)
Zmax = max(Z.max(), abs(Z.min()))

norm = colors.SymLogNorm(
    linthresh=1e-3,   # size of linear region around 0
    linscale=1,
    vmin=-Zmax,
    vmax=Zmax
)
surf= ax.contour(X, Y, Z, levels=400, norm=norm, cmap="RdBu_r")
#surf = ax.streamplot(X, Y, Z_x, Z_y, density=5.)
plt.show()
