import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# load data
data = np.loadtxt("pole_evolution.csv", delimiter=",")

t = data[:, 0]
zdata = data[:, 1:]
n_poles = zdata.shape[1] // 2
n_steps = len(t)

# prepare figure
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect("equal")
ax.set_xlabel("Re(z)")
ax.set_ylabel("Im(z)")
ax.set_title("Pole trajectories")

# plot limits (fixed, important for smooth animation)
ax.set_xlim(np.min(zdata[:, ::2]), np.max(zdata[:, ::2]))
ax.set_ylim(np.min(zdata[:, 1::2]), np.max(zdata[:, 1::2]))

# one line per pole
lines = []
for _ in range(n_poles):
    line, = ax.plot([], [], lw=1)
    lines.append(line)

def update(frame):
    for k, line in enumerate(lines):
        re = zdata[:frame, 2*k]
        im = zdata[:frame, 2*k + 1]
        line.set_data(re, im)
    ax.set_title(f"time = {t[frame]:.2f}")
    return lines

ani = FuncAnimation(
    fig,
    update,
    frames=n_steps,
    interval=30,   # ms between frames
    blit=False,
    repeat=False
)

plt.show()

