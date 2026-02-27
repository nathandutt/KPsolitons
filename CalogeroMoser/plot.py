#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load data
A = np.load("Output/cpoints.npy")   # shape (T, N, 2)
T, N, _ = A.shape

# --- Compute global limits ignoring NaNs ---
xmin = np.nanmin(A[:, :, 0])
xmax = np.nanmax(A[:, :, 0])
ymin = np.nanmin(A[:, :, 1])
ymax = np.nanmax(A[:, :, 1])
ymin -=1
ymax+=1

# --- Add padding (5%) ---
pad_frac = 0.05

xrange = xmax - xmin
yrange = ymax - ymin

xmin -= pad_frac * xrange
xmax += pad_frac * xrange
ymin -= pad_frac * yrange
ymax += pad_frac * yrange

fig, ax = plt.subplots()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# Optional: equal aspect ratio (usually prettier for particle systems)
ax.set_aspect('equal', adjustable='box')

# Initial frame (remove NaNs)
pts0 = A[0]
pts0 = pts0[~np.isnan(pts0).any(axis=1)]

scat = ax.scatter(pts0[:, 0], pts0[:, 1], s=10, alpha=0.7)

def update(frame):
    pts = A[frame]
    pts = pts[~np.isnan(pts).any(axis=1)]
    scat.set_offsets(pts)
    return scat,

ani = FuncAnimation(
    fig,
    update,
    frames=T,
    interval=2000//T,
    blit=True
)

plt.show()
