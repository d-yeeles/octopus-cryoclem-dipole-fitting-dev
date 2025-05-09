import numpy as np
from MLEwT_fixed_2 import dipdistr
from math import ceil
import os
import sys
import time
import re
import datetime
from vectorized_mortensen import FullyVectorizedDipolePSFGenerator

class DipolePSFGenerator:
    def __init__(self, image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor=3, verbose=False):
        self.image_size = (int(image_size[0]), int(image_size[1]))
        self.subpixel_factor = subpixel_factor
        self.pixel_size = pixel_size
        self.wavelength = wavelength
        self.n_sample = n_sample
        self.n_objective = n_objective
        self.magnification = magnification
        self.NA = NA
        self.norm_file = norm_file
        self.verbose = verbose

        # Sort out subpixel averaging bits
        self.oversampled_pixel_size = self.pixel_size/self.subpixel_factor
        self.oversampled_image_size = (self.image_size[0]*self.subpixel_factor, self.image_size[1]*self.subpixel_factor)
        
        # Load normalization file
        self.norm_data = np.load(norm_file)
        
        # Initialize dipole distribution model
        self.DD = dipdistr(wavelength, n_objective, n_sample, magnification, NA, norm_file)
    
    def __call__(self, phi, theta, x_pos, y_pos, n_photons):

        # Create position vector
        # x_pos and y_pos are in nm
        posvec = np.arange(-(self.oversampled_image_size[0]-1)/2, self.oversampled_image_size[0]/2) * self.oversampled_pixel_size
        dipole_psf = np.zeros(self.oversampled_image_size)

        # Generate PSF
#        start_time_1 = time.time()
        for i in range(self.oversampled_image_size[0]):
            for j in range(self.oversampled_image_size[1]):
         
#                start_time_2 = time.time()
                dipole_psf[j, i] = self.DD.PSF_approx(posvec[i] - x_pos, 
                                                      -posvec[j] - y_pos,
                                                      phi, theta)
#                end_time_2 = time.time()
#                elapsed_time_2 = end_time_2 - start_time_2
#                print(f"Subpixel: {elapsed_time_2:.6f} seconds")
#        end_time_1 = time.time()
#        elapsed_time_1 = end_time_1 - start_time_1
#        print(f"Whole image: {elapsed_time_1:.6f}")
       
#        dipole_psf_upsampled = dipole_psf.copy()
#        dipole_psf_upsampled = dipole_psf_upsampled / dipole_psf_upsampled.sum() * n_photons

        # Return back to the real coarser pixel size after (subpixel averaging)
        dipole_psf = dipole_psf.reshape(self.image_size[0], self.subpixel_factor, 
                                     self.image_size[1], self.subpixel_factor)
        dipole_psf = np.mean(dipole_psf, axis=(1, 3))
#        dipole_psf = dipole_psf / dipole_psf.sum() * n_photons

#        print(f"Upsampled sum: {dipole_psf_upsampled.sum()}, Downsampled sum: {dipole_psf.sum()}")

#        import tifffile
#        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
#        filename = f"dipole_psf_{timestamp}.tiff"
#        # Normalize to 0-65535 for 16-bit TIFF
#        if dipole_psf.min() != dipole_psf.max():
#            normalized = ((dipole_psf - dipole_psf.min()) / (dipole_psf.max() - dipole_psf.min()) * 65535).astype(np.uint16)
#        else:
#            normalized = np.zeros_like(dipole_psf, dtype=np.uint16)
#        tifffile.imwrite(filename, normalized)
#        print(f"Saved dipole PSF to {filename}")

        return dipole_psf



def run_simulator(x, y, theta, phi, image_size_px, pixel_size, wavelength, n_objective, n_sample, NA, n_photons):
    """
    Runs the Mortensen fit for given phi and theta, returning the results and ground truth.
    """

    # Define parameters
    image_size = (image_size_px, image_size_px)  
    magnification = 215 
    norm_file = "/home/tfq96423/dipolenorm.npy"

#    # Define dipole ground_truth position in nm
#    x_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2)  
#    y_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2) 
#    
    # factor by which to oversample
    subpixel_factor = 1

    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor)
    dipole_psf = psf_generator(phi, theta, x, y, n_photons)

#    # Vectorised version
#    psf_generator = FullyVectorizedDipolePSFGenerator(
#        image_size, pixel_size, wavelength, n_objective, n_sample, 
#        magnification, NA, norm_file, subpixel_factor
#    )
#    start_time_1 = time.time()
#    dipole_psf = psf_generator(phi, theta, x, y, n_photons)
#    end_time_1 = time.time()
#    elapsed_time_1 = end_time_1 - start_time_1
#    print(f"Whole image: {elapsed_time_1:.6f}")

    # Generate dipole PSF
    # this is the slow part
#    start_time = time.time()
#    dipole_psf = psf_generator(phi, theta, x, y, n_photons)
#    end_time = time.time()
#    elapsed_time = end_time - start_time
#    print(f"slow bit: {elapsed_time:.2f}")

#    dipole_psf_noisy = np.random.poisson(dipole_psf)

    return dipole_psf#_noisy

