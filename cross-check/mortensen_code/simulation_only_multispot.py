import numpy as np
from MLEwT_fixed import dipdistr, MLEwT
from math import ceil
import os
import sys
import time
import re
import datetime
#import matplotlib.pyplot as plt
import tifffile
from save_results import save_results_to_py

def parse_parameters_file(file_path):
    import re
    settings = {}
    # Open the file and read the lines
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Define regex patterns for different types of lines
    number_pattern = re.compile(r'([\w\.]+)\s*=\s*(-?\d+\.?\d*)')  # Updated to include periods in parameter names
    array_pattern = re.compile(r'([\w\.]+)\s*=\s*\[(.*)\]')       # Updated to include periods in parameter names
    # Loop through each line to process it
    for line in lines:
        line = line.strip()
        # Skip comment lines
        if line.startswith('%') or not line:
            continue
        # Match number assignments
        match = number_pattern.match(line)
        if match:
            var_name, var_value = match.groups()
            settings[var_name] = float(var_value) if '.' in var_value else int(var_value)
            continue
        # Match array assignments
        match = array_pattern.match(line)
        if match:
            var_name, array_values = match.groups()
            # Convert the array string into a list of floats, handling empty entries
            values = []
            for val in array_values.split(','):
                val = val.strip()
                if val:  # Only process non-empty values
                    values.append(float(val))
            settings[var_name] = values
            continue
    return settings

def save_as_tif(psf_total_image, output_path):
    # Ensure the image is uint32
    psf_total_image = psf_total_image.astype(np.uint32)
    
    # Save as TIF with LZW compression
    with tifffile.TiffWriter(output_path, bigtiff=True) as t:
        t.write(
            psf_total_image, 
            photometric='minisblack',  # Grayscale
            planarconfig='contig',  # Add this to match MATLAB's default
            dtype=np.uint32,  # 32-bit integer
            rowsperstrip=16,  # Strip length
            compression='lzw',  # Lossless compression
            software='Python'  # Software metadata
        )


class DipolePSFGenerator:
    def __init__(self, image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor=3, verbose=False):
        self.image_size = image_size
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
        for i in range(self.oversampled_image_size[0]):
            for j in range(self.oversampled_image_size[1]):
                dipole_psf[j, i] = self.DD.PSF_approx(posvec[i] - x_pos, 
                                                      posvec[j] - y_pos,
                                                      phi, theta, 
                                                      )
        
#        dipole_psf_upsampled = dipole_psf.copy()
#        dipole_psf_upsampled = dipole_psf_upsampled / dipole_psf_upsampled.sum() * n_photons


        # Return back to the real coarser pixel size after (subpixel averaging)
        dipole_psf = dipole_psf.reshape(self.image_size[0], self.subpixel_factor, 
                                     self.image_size[1], self.subpixel_factor)
        dipole_psf = np.mean(dipole_psf, axis=(1, 3))

        dipole_psf = dipole_psf / dipole_psf.sum() * n_photons

#        print(f"Upsampled sum: {dipole_psf_upsampled.sum()}, Downsampled sum: {dipole_psf.sum()}")

        return dipole_psf

def run_simulator(phi, theta, grid_position=None, grid_dims=None, image_size=None):
    """
    Runs the Mortensen fit for given phi and theta, returning the results and ground truth.
    """
    # Define parameters
    if image_size is None:
        image_size = (19, 19)
    pixel_size = 52.0
    wavelength = 500
    n_objective = 2.17  
    n_sample = 1.31  
    magnification = 215 
    NA = 2.17  
    norm_file = "/home/tfq96423/dipolenorm.npy"
    n_photons = 2000
    number_of_spots = 1
    background = 0
    objectiveFocalLength = 770
    nDiscretizationBFP = 129

    # Define dipole ground_truth position in nm (random)
    x_pos_nm_local = np.random.uniform(-image_size[1]/2, image_size[1]/2)  
    y_pos_nm_local = np.random.uniform(-image_size[0]/2, image_size[0]/2) 

    
    # Calculate global position if grid information is provided
    if grid_position is not None and grid_dims is not None:
        row, col = grid_position
        grid_rows, grid_cols = grid_dims
        
        # Calculate full grid size in nm
        grid_width_nm = grid_cols * image_size[1] * pixel_size
        grid_height_nm = grid_rows * image_size[0] * pixel_size
        
        # Calculate patch centre position in the grid (in nm, relative to grid centre)
        patch_centre_x_nm = (col + 0.5) * image_size[1] * pixel_size - grid_width_nm/2
        patch_centre_y_nm = (row + 0.5) * image_size[0] * pixel_size - grid_height_nm/2
        
        # Calculate global position by adding local offset to patch centre
        x_pos_nm_global = patch_centre_x_nm + x_pos_nm_local
        y_pos_nm_global = patch_centre_y_nm + y_pos_nm_local
    else:
        # If no grid information, just use local positions
        x_pos_nm_global = x_pos_nm_local
        y_pos_nm_global = y_pos_nm_local


    # factor by which to oversample
    subpixel_factor = 1

    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor)

    # Generate dipole PSF
    dipole_psf = psf_generator(phi, theta, x_pos_nm_local, y_pos_nm_local, n_photons)
    dipole_psf_noisy = np.random.poisson(dipole_psf)

    parameters = [number_of_spots, pixel_size, wavelength, NA, objectiveFocalLength, [n_sample, n_objective, n_objective], nDiscretizationBFP, background, n_photons, x_pos_nm_global, y_pos_nm_global, theta, phi]

    # Return results instead of printing
    return dipole_psf_noisy, parameters

# Prevent script from running when imported
if __name__ == "__main__":

    number_of_dipoles = 1

    runs = range(1)
    thetas = [67.5*np.pi/180]#np.arange(0, np.pi/2 + 0.001, 22.5*np.pi/180)#np.arange(0, np.pi/2 + 0.001, 22.5*np.pi/180)
    phis = np.arange(0, 2*np.pi - 0.001, 0.25*np.pi/180)

    for run in runs:
      for theta in thetas:
        for phi in phis:

            theta_deg = theta*180/np.pi
            phi_deg = phi*180/np.pi

            # Define output paths
#            image_output_path = f'../sims_1spot_loads/sim_theta{ceil(theta_deg):03}_phi{ceil(phi_deg):03}_run{int(run+1)}.tif'
            image_output_path = f'../sims_1spot_loads/sim_theta{theta_deg:.2f}_phi{phi_deg:.2f}_run{int(run+1)}.tif'
#            params_output_path = f"../sims_1spot_loads/params_theta{ceil(theta_deg):03}_phi{ceil(phi_deg):03}_run{int(run+1)}.m"
            params_output_path = f'../sims_1spot_loads/params_theta{theta_deg:.2f}_phi{phi_deg:.2f}_run{int(run+1)}.m'

            # Initialise arrays to collect parameters for all dipoles
            positionX_nm_list = []
            positionY_nm_list = []
            angleInclination_list = []
            angleAzimuth_list = []
                
            # Initialise the total image
            #psf_total_image = None
            # Calculate grid dimensions
            grid_cols = ceil(np.sqrt(number_of_dipoles))
            grid_rows = ceil(number_of_dipoles / grid_cols)

            # Calculate grid dimensions
            grid_cols = ceil(np.sqrt(number_of_dipoles))
            grid_rows = ceil(number_of_dipoles / grid_cols)
            grid_dims = (grid_rows, grid_cols)

            # Get image dimensions from the first run
            sample_image, _ = run_simulator(phi, theta)
            image_height, image_width = sample_image.shape

            # Create an empty grid image
            psf_total_image = np.zeros((grid_rows * image_height, grid_cols * image_width), dtype=sample_image.dtype)



            # Loop over each dipole in the frame
            for dipole in range(number_of_dipoles):
                # Calculate position in grid
                row = dipole // grid_cols
                col = dipole % grid_cols
                grid_position = (row, col)

                # Generate the PSF image and parameters with global coordinates
                psf_image, parameters = run_simulator(
                    phi, theta, 
                    grid_position=grid_position, 
                    grid_dims=grid_dims,
                    image_size=(image_height, image_width)
                )

                # Collect parameters for this dipole
                positionX_nm_list.append(parameters[9])  # Now contains global X coordinate
                positionY_nm_list.append(parameters[10])  # Now contains global Y coordinate
                angleInclination_list.append(parameters[11])
                angleAzimuth_list.append(parameters[12])
                          
                # Place the image in the grid
                y_start = row * image_height
                y_end = (row + 1) * image_height
                x_start = col * image_width
                x_end = (col + 1) * image_width
                psf_total_image[y_start:y_end, x_start:x_end] = psf_image
        
            # Save total image
            save_as_tif(psf_total_image, image_output_path)
     
            # Write data to file
            with open(params_output_path, 'w') as file:
                file.write(f"% ground truth for {image_output_path}\n")
                file.write("% settings\n")
                file.write(f"number_of_spots = {parameters[0]}\n")
                file.write(f"pixel_size_nm = {parameters[1]}\n")
                file.write(f"image_size_nm = {psf_total_image.shape[0] * parameters[1]}\n")
                file.write(f"image_size_px = {psf_total_image.shape[0]}\n")
                file.write(f"wavelength = {parameters[2]}\n")
                file.write(f"par.objectiveNA = {parameters[3]}\n")
                file.write(f"objectiveFocalLength = {parameters[4]}\n")
                file.write(f"par.refractiveIndices = {parameters[5]}\n")
                file.write(f"par.nDiscretizationBFP = {parameters[6]}\n")
                file.write(f"par.backgroundNoise = {parameters[7]}\n")
                file.write(f"par.nPhotons = {parameters[8]}\n")
                file.write(f"positionX_nm_array = [{', '.join(str(x) for x in positionX_nm_list)}]\n")
                file.write(f"positionY_nm_array = [{', '.join(str(y) for y in positionY_nm_list)}]\n")
                file.write(f"angleInclination_array = [{', '.join(str(theta) for theta in angleInclination_list)}]\n")
                file.write(f"angleAzimuth_array = [{', '.join(str(phi) for phi in angleAzimuth_list)}]\n")
            file.close()

