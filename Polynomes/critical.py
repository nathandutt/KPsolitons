import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ======================
# USER PARAMETERS
# ======================
filename = "coefficients.txt"

l = 2.0
x_min, x_max = -10 * l, 10 * l
y_min, y_max = -10 * l, 10 * l
nx, ny = 2000, 2000   # slightly reduced for safety

eps = 1e-8
max_points = 3000

# ======================
# POLYNOMIAL UTILITIES
# ======================
def reconstruct_polynomial(coeffs, X, Y):
    Z = np.zeros_like(X)
    idx = 0
    d = 0
    while idx < len(coeffs):
        for n in range(d + 1):
            if idx >= len(coeffs):
                break
            Z += coeffs[idx] * (X ** n) * (Y ** (d - n))
            idx += 1
        d += 1
    return Z


def reconstruct_gradients(coeffs, X, Y):
    dPx = np.zeros_like(X)
    dPy = np.zeros_like(X)

    idx = 0
    d = 0
    while idx < len(coeffs):
        for n in range(d + 1):
            if idx >= len(coeffs):
                break
            c = coeffs[idx]
            px = n
            py = d - n

            if px > 0:
                dPx += c * px * (X ** (px - 1)) * (Y ** py)
            if py > 0:
                dPy += c * py * (X ** px) * (Y ** (py - 1))

            idx += 1
        d += 1

    return dPx, dPy


# ======================
# READ FILE
# ======================
times = []
P_coeffs = []
Q_coeffs = []

with open(filename, "r") as f:
    lines = f.readlines()

for i in range(0, len(lines), 2):
    t1, *p = map(float, lines[i].split())
    t2, *q = map(float, lines[i + 1].split())
    assert t1 == t2
    times.append(t1)
    P_coeffs.append(p)
    Q_coeffs.append(q)

times = np.array(times)

# ======================
# GRID
# ======================
x = np.linspace(x_min, x_max, nx)
y = np.linspace(y_min, y_max, ny)
X, Y = np.meshgrid(x, y)

# ======================
# LOCAL MINIMA (NO THRESHOLD)
# ======================
def local_minima(G):
    return (
        (G < np.roll(G,  1, axis=0)) &
        (G < np.roll(G, -1, axis=0)) &
        (G < np.roll(G,  1, axis=1)) &
        (G < np.roll(G, -1, axis=1)) &
        (G < np.roll(G, (1, 1), axis=(0, 1))) &
        (G < np.roll(G, (1,-1), axis=(0, 1))) &
        (G < np.roll(G, (-1,1), axis=(0, 1))) &
        (G < np.roll(G, (-1,-1), axis=(0, 1)))
    )

# ======================
# FIGURE
# ======================
fig, ax = plt.subplots()
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_xlabel("x")
ax.set_ylabel("y")

sc = ax.scatter([], [], s=4, c="black")
title = ax.set_title("")

# ======================
# UPDATE
# ======================
def update(frame):
    P = reconstruct_polynomial(P_coeffs[frame], X, Y)
    Q = reconstruct_polynomial(Q_coeffs[frame], X, Y)

    Px, Py = reconstruct_gradients(P_coeffs[frame], X, Y)
    Qx, Qy = reconstruct_gradients(Q_coeffs[frame], X, Y)

    Gx = P * Qx - Q * Px
    Gy = P * Qy - Q * Py
    G = np.sqrt(Gx**2 + Gy**2)

    G[np.abs(P) < eps] = np.inf

    mask = local_minima(G)

    xs = X[mask]
    ys = Y[mask]

    print(f"frame {frame}: {xs.size} minima")

    if xs.size > max_points:
        idx = np.random.choice(xs.size, max_points, replace=False)
        xs = xs[idx]
        ys = ys[idx]

    sc.set_offsets(np.c_[xs, ys])
    title.set_text(f"Local minima of |∇(Q/P)|   t = {times[frame]:.3f}")

    return sc, title


ani = FuncAnimation(fig, update, frames=len(times), interval=80)
plt.show()

