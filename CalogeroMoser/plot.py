import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Directory containing .npy files
directory = "Output"

# Collect (time, filepath) pairs
time_file_pairs = []

for filename in os.listdir(directory):
    if filename.endswith(".npy"):
        name_without_ext = os.path.splitext(filename)[0]
        try:
            t = float(name_without_ext)
            filepath = os.path.join(directory, filename)
            time_file_pairs.append((t, filepath))
        except ValueError:
            # Skip files that don't match float naming
            pass

# Sort by time
time_file_pairs.sort(key=lambda x: x[0])

# Extract sorted times and filepaths
times = [pair[0] for pair in time_file_pairs]
files = [pair[1] for pair in time_file_pairs]

# Preload all data
data_list = [np.load(f) for f in files]

# Determine global axis limits for consistent animation scaling
all_points = np.vstack(data_list)
xmin, ymin = np.min(all_points, axis=0)
xmax, ymax = np.max(all_points, axis=0)

# Create figure
fig, ax = plt.subplots()
scat = ax.scatter([], [])

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
title = ax.set_title("")

def init():
    scat.set_offsets(np.empty((0, 2)))
    title.set_text("")
    return scat, title

def update(frame):
    points = data_list[frame]
    scat.set_offsets(points)
    title.set_text(f"Time = {times[frame]:.4f}")
    return scat, title

ani = FuncAnimation(
    fig,
    update,
    frames=len(data_list),
    init_func=init,
    interval=50,  # milliseconds between frames
    blit=True
)
ani.save("aniamtion.gif", writer="ffmpeg", fps=25)
plt.show()
