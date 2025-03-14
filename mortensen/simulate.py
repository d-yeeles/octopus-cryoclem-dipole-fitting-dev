import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

# Import the original implementation
from diPOLE_python3 import dipdistr as dipdistr_original

# Import our fixed HintererPSF
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diPOLE_python3_hinterer_auto import dipdistr as dipdistr_hinterer

# Set common parameters
peak_emission_wavelength = 580.0  # nm
ref_ind_immersion = 1.52
ref_ind_buffer = 1.33
numerical_aperture = 1.49
magnification = 250.0
pixel_size_nm = 51.2/10


def plot_original_psf_grid():
    """Generate and display raw PSF using the original implementation for theta from 0 to 90 degrees"""
    # Create the dipole distribution object
    dd = dipdistr_original(peak_emission_wavelength, ref_ind_immersion, ref_ind_buffer, magnification,
                           numerical_aperture)

    # Common parameters
    phi = 0.0  # Azimuthal angle (rad)
    deltaz = 0.0  # Defocus

    # Generate grid
    image_size_px = 170
    image_size_nm = image_size_px * pixel_size_nm
    pixels = 170

    x = np.linspace(-image_size_nm / 2, image_size_nm / 2, pixels)
    y = np.linspace(-image_size_nm / 2, image_size_nm / 2, pixels)
    X, Y = np.meshgrid(x, y)

    # Create figure for 4x2 grid
    fig, axes = plt.subplots(1, 3, figsize=(16, 10))
    axes = axes.flatten()

    # Generate 8 different theta values from 0 to 90 degrees
    theta_values = np.linspace(0, np.pi / 2, 3)

    # Calculate and plot PSF for each theta value
    for i, theta in enumerate(theta_values):
        # Calculate PSF
        Z = np.zeros_like(X)
        for ii in range(pixels):
            for jj in range(pixels):
                Z[ii, jj] = dd.PSF_approx(X[ii, jj], Y[ii, jj], phi, theta, deltaz)

        # Plot on corresponding subplot
        im = axes[i].imshow(Z, cmap='gray', extent=[-image_size_nm/2, image_size_nm/2, -image_size_nm/2, image_size_nm/2])
        axes[i].set_title(f'θ={np.degrees(theta):.1f}°')
        axes[i].set_xlabel('x (nm)')
        axes[i].set_ylabel('y (nm)')
        # Add grid for better readability
        axes[i].grid(False)

    plt.tight_layout()
    plt.show()
    plt.close()


def plot_hinterer_psf_grid():
    """Generate and display raw PSF using the Hinterer implementation for theta from 0 to 90 degrees"""
    # Create the dipole distribution object
    dd = dipdistr_hinterer(peak_emission_wavelength, ref_ind_immersion, ref_ind_buffer, magnification,
                           numerical_aperture)

    # Common parameters
    phi = 0.0  # Azimuthal angle (rad)
    deltaz = 0.0  # Defocus

    # Generate grid
    image_size_px = 170
    image_size_nm = image_size_px * pixel_size_nm
    pixels = 170

    x = np.linspace(-image_size_nm / 2, image_size_nm / 2, pixels)
    y = np.linspace(-image_size_nm / 2, image_size_nm / 2, pixels)
    X, Y = np.meshgrid(x, y)

    # Create figure for 4x2 grid
    fig, axes = plt.subplots(2, 4, figsize=(16, 10))
    axes = axes.flatten()

    # Generate 8 different theta values from 0 to 90 degrees
    theta_values = np.linspace(0, np.pi / 2, 8)

    # Calculate and plot PSF for each theta value
    for i, theta in enumerate(theta_values):
        # Calculate PSF
        Z = np.zeros_like(X)
        for ii in range(pixels):
            for jj in range(pixels):
                Z[ii, jj] = dd.PSF_approx(X[ii, jj], Y[ii, jj], phi, theta, deltaz)
                # Z = dd.PSF_approx(X[ii, jj], Y[ii, jj], phi, theta, deltaz)
        # Plot on corresponding subplot with nm coordinates
        im = axes[i].imshow(Z, cmap='gray', extent=[-image_size_nm/2, image_size_nm/2, -image_size_nm/2, image_size_nm/2])
        axes[i].set_title(f'θ={np.degrees(theta):.1f}°')
        axes[i].set_xlabel('x (nm)')
        axes[i].set_ylabel('y (nm)')
        # Add grid for better readability
        axes[i].grid(False)

    plt.tight_layout()
    plt.show()
    plt.close()


if __name__ == "__main__":
    # Generate grid plots for both implementations
    print("Generating original PSF grid...")
    # plot_original_psf_grid()

    print("Generating Hinterer PSF grid...")
    plot_hinterer_psf_grid()