import numpy as np
import matplotlib.pyplot as plt

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
x_s = np.linspace(-x_range, x_range, 1000)

def field(x, Z_s):
    return np.real(np.sum(1 / (x - Z_s)))

# Precompute phi
phi = np.zeros((n_steps, x_s.size))
for i in range(n_steps):
    for j in range(x_s.size):
        phi[i, j] = field(x_s[j], all_Z[i])

# ---- Animation state ----
frame = 0
playing = False

fig, ax = plt.subplots(figsize=(6, 6))
(line,) = ax.plot(x_s, phi[frame])
ax.set_ylim(-10, 10)
ax.set_title("Hold SPACE to play")

def update():
    global frame
    if playing:
        frame = (frame + 1) % n_steps
        line.set_ydata(phi[frame])
        fig.canvas.draw_idle()

def on_key_press(event):
    global playing
    if event.key == " ":
        playing = True

def on_key_release(event):
    global playing
    if event.key == " ":
        playing = False

fig.canvas.mpl_connect("key_press_event", on_key_press)
fig.canvas.mpl_connect("key_release_event", on_key_release)

timer = fig.canvas.new_timer(interval=30)
timer.add_callback(update)
timer.start()

plt.show()

