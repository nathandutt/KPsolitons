import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import SymLogNorm
from matplotlib.cm import ScalarMappable

# ======================
# USER PARAMETERS
# ======================
filename = "coefficients.txt"
l = 1.
x_min, x_max = l*-14.0, l*4.0
y_min, y_max = -l*7.0, l*7.0
nx, ny = 600, 600

levels = 200
eps = 1e-9
cmap = "coolwarm"

linthresh = 1e-8
linscale = 1.0

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


def symlog_levels(vmax, linthresh, n):
    pos = np.logspace(
        np.log10(linthresh),
        np.log10(vmax),
        n // 2
    )
    return np.concatenate((-pos[::-1], [0.0], pos))


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
# INITIAL DATA SCALING
# ======================
P0 = reconstruct_polynomial(P_coeffs[0], X, Y)
Q0 = reconstruct_polynomial(Q_coeffs[0], X, Y)
Z0 = np.where(np.abs(P0) > eps, Q0 / P0, np.nan)

Z0_finite = Z0[np.isfinite(Z0)]
if Z0_finite.size == 0:
    raise RuntimeError("First frame has no finite values.")

vmax = np.nanpercentile(np.abs(Z0_finite), 99)
vmax = max(vmax, linthresh * 10)
vmin = -vmax

norm = SymLogNorm(
    linthresh=linthresh,
    linscale=linscale,
    vmin=vmin,
    vmax=vmax
)

levels_arr = symlog_levels(vmax, linthresh, levels)

# ======================
# FIGURE
# ======================
fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("y")

# Initial contours (may be empty — that's OK)
contour = ax.contour(
    X, Y, Z0,
    levels=levels_arr,
    cmap=cmap,
    norm=norm
)

# SAFE COLORBAR (independent of contours)
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # required by matplotlib
cb = plt.colorbar(sm, ax=ax)
cb.set_label("Q / P")

title = ax.set_title(f"t = {times[0]:.3f}")

# ======================
# ANIMATION FUNCTION
# ======================
def update(frame):
    global contour

    for c in contour.collections:
        c.remove()

    P = reconstruct_polynomial(P_coeffs[frame], X, Y)
    Q = reconstruct_polynomial(Q_coeffs[frame], X, Y)
    Z = np.where(np.abs(P) > eps, Q / P, np.nan)

    contour = ax.contour(
        X, Y, Z,
        levels=levels_arr,
        cmap=cmap,
        norm=norm
    )

    title.set_text(f"t = {times[frame]:.3f}")
    return contour.collections + [title]


ani = FuncAnimation(
    fig,
    update,
    frames=len(times),
    interval=50,
    blit=False,
    repeat=False,
)

# ======================
# SAVE OR SHOW
# ======================
# ani.save("QP_symlog_contours.mp4", writer="ffmpeg", fps=10, dpi=150)
plt.show()

