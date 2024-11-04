import numpy as np
import matplotlib.pyplot as plt

def gaussian_kernel(size, sigma=1):
    """Generate a Gaussian kernel."""
    kernel = np.fromfunction(
        lambda x, y: (1 / (2 * np.pi * sigma ** 2)) *
                     np.exp(-((x - (size - 1) / 2) ** 2 + (y - (size - 1) / 2) ** 2) / (2 * sigma ** 2)),
        (size, size)
    )
    return kernel / np.sum(kernel)  # Normalize the kernel

# Generate a 5x5 Gaussian kernel with sigma = 1
size = 5

sigma = 1
kernel = gaussian_kernel(size, sigma)

# Display the kernel
plt.imshow(kernel, cmap='gray')
plt.title('Gaussian Filter (σ = 1.0)')
plt.colorbar()
plt.show()