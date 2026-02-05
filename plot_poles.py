import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ----------------------
# load data
# ----------------------
data = np.loadtxt("pole_evolution.csv", delimiter=",")

t = data[:, 0]
zdata = data[:, 1:]
n_poles = zdata.shape[1] // 2
n_steps = len(t)
timestep = 0.2
# ----------------------
# prepare figure
# ----------------------
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect("equal")
ax.set_xlabel("Re(z)")
ax.set_ylabel("Im(z)")
ax.set_title("Pole trajectories (click to play/pause)")
y_range = 2*(n_steps*timestep)
ax.set_xlim(-y_range, y_range)
ax.set_ylim(-y_range, y_range)

# one line per pole
lines = []
for _ in range(n_poles):
    line, = ax.plot([], [], marker ='o')
    lines.append(line)

# ----------------------
# update function
# ----------------------
def update(frame):
    for k, line in enumerate(lines):
        re = zdata[:frame, 2 * k]
        im = zdata[:frame, 2 * k + 1]
        line.set_data(re, im)

    ax.set_title(f"y = {t[frame]:.2f}")
    return lines

# ----------------------
# create animation
# ----------------------
ani = FuncAnimation(
    fig,
    update,
    frames=n_steps,
    interval=5,
    blit=False,
    repeat=True
)

# start paused
ani.event_source.stop()

# ----------------------
# mouse click handler
# ----------------------
playing = False

def on_click(event):
    global playing
    if event.button != 1:  # left mouse button only
        return

    if playing:
        ani.event_source.stop()
        playing = False
        ax.set_title("Pole trajectories (paused)")
    else:
        ani.event_source.start()
        playing = True

fig.canvas.mpl_connect("button_press_event", on_click)

# ----------------------
# show
# ----------------------
plt.show()

