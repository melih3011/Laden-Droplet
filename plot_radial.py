#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sys import argv

def centeredDistanceMatrix(n):
    """
    Generates a distance matrix with distances centered around the middle.
    
    Parameters:
    n (int): Size of the matrix (should be odd).
    
    Returns:
    np.ndarray: Distance matrix.
    """
    if n % 2 == 0:
        raise ValueError("n should be odd.")
    x, y = np.meshgrid(np.arange(n), np.arange(n))
    return np.sqrt((x - (n // 2))**2 + (y - (n // 2))**2)

def arbitraryfunction(d, y):
    """
    Interpolates values based on distance matrix and provided data.
    
    Parameters:
    d (np.ndarray): Distance matrix.
    y (np.ndarray): Values to interpolate.
    
    Returns:
    np.ndarray: Interpolated values reshaped to the shape of d.
    """
    n = len(y)
    x = np.arange(n)
    f = interp1d(x, y, bounds_error=False, fill_value="none")
    return f(d.flat).reshape(d.shape)

# Load data from file
data = np.loadtxt(argv[1])

# Ensure data has at least 2 columns
if data.shape[1] < 2:
    raise ValueError("Data file must have at least two columns.")

# Extract the columns
y1 = data[:, 1]  # First column
y2 = data[:, 2]  # Second column

# Set n to the length of y (assumes y1 and y2 have the same length)
n = len(y1)

# Optional: Handle even n by making it odd
if n % 2 == 0:
    print(f"Warning: Length of data ({n}) is even. Adjusting to an odd value.")
    n -= 1
    y1 = y1[:n]
    y2 = y2[:n]

# Generate distance matrix and interpolated values
d = centeredDistanceMatrix(n)
f1 = arbitraryfunction(d, y1)
f2 = arbitraryfunction(d, y2)

# Displaying the 2D interpolated matrices
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Display the 2D matrix for the first column
cax1 = axs[0].imshow(f1.T, origin='lower', interpolation='antialiased', cmap='viridis')
fig.colorbar(cax1, ax=axs[0])
axs[0].set_title('Liquid')
axs[0].set_xticks([])
axs[0].set_yticks([])

# Display the 2D matrix for the second column
cax2 = axs[1].imshow(f2.T, origin='lower', interpolation='antialiased', cmap='viridis')
fig.colorbar(cax2, ax=axs[1])
axs[1].set_title('Nano')
axs[1].set_xticks([])
axs[1].set_yticks([])

plt.tight_layout()
plt.savefig('output.png', dpi=300)
plt.show()
