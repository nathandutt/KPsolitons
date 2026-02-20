import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Directory containing .npy files
output_dir = "Output"
npy_files = sorted(glob.glob(os.path.join(output_dir, "*.npy")))

plt.figure(figsize=(10, 6))

for npy_file in npy_files:
    # Load data: shape (Ny, Ns, 2)
    data = np.load(npy_file)
    Ny, Ns, _ = data.shape

    # Plot each of the Ns lines (x,y) for all Ny parameter values
    for i in range(Ns):
        x = data[:, i, 0]  # x-coordinates for line i
        y = data[:, i, 1]  # y-coordinates for line i
        plt.plot(x, y, '-', linewidth=1, label=f"{os.path.basename(npy_file)} (line {i})" if i == 0 else "")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Lines from all .npy files (Ns lines per file, traced over Ny parameter values)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()

