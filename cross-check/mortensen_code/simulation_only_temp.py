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





    def mortensen_fit(self, dipole_psf):
 
        # Define parameters
        #deltapix = 9
        
        # Generate the initial mux, muy in pixels around the centre of the image (9, 9)
        mux_pix = np.random.uniform(self.image_size[1] / 2 - 0.5, self.image_size[1]/ 2 + 0.5)
        muy_pix = np.random.uniform((self.image_size[0] - 1)/ 2 - 0.5, (self.image_size[0] - 1) / 2 + 0.5)
        #print(f"the initial vals of mux muy: {mux_pix, muy_pix}")

        # Convert to coordinates in nm
        mux_nm = np.random.uniform(0 - self.pixel_size/2, 0 + self.pixel_size/2)    
        muy_nm = np.random.uniform(0 - self.pixel_size/2, 0 + self.pixel_size/2)    
        
#        init_theta = np.random.uniform(0, np.pi/2)
#        init_phi = np.random.uniform(0, 2 * np.pi)
        init_new1 = np.random.uniform(-1, 1)
        init_new2 = np.random.uniform(-1, 1)
        init_new3 = np.random.uniform(0, 1)

        init_photons = np.sum(dipole_psf)
        
#        initvals = np.array([init_phi, init_theta, mux_nm, muy_nm, init_photons])
        initvals = np.array([init_new1, init_new2, init_new3, mux_nm, muy_nm, init_photons])
        
        #initpix = (self.image_size[0] // 2, self.image_size[1] // 2) # (ypix, xpix)
        #print(f"initial values for the centre pixel (ypixel,xpixel): {initpix}")
        
        # initvals = the initial values to start the estimate       
        track = MLEwT(initvals, self)   #initpix, deltapix, self)
        result = track.Estimate(dipole_psf)

        return result


def run_simulator(phi, theta):
    """
    Runs the Mortensen fit for given phi and theta, returning the results and ground truth.
    """
    # Define parameters
    image_size = (155, 155)  
    pixel_size = 52.0/8.0
    wavelength = 500
    n_objective = 2.17  
    n_sample = 1.31  
    magnification = 215 
    NA = 2.17  
    norm_file = "/home/tfq96423/dipolenorm.npy"
    n_photons = 1e10
    number_of_spots = 1
    background = 0
    objectiveFocalLength = 770
    nDiscretizationBFP = 129

    # Define dipole ground_truth position in nm
    x_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2)  
    y_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2) 
    
    # factor by which to oversample
    subpixel_factor = 1

    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor)

    # Generate dipole PSF
    dipole_psf = psf_generator(phi, theta, x_pos_nm, y_pos_nm, n_photons)
    dipole_psf_noisy = np.random.poisson(dipole_psf)

    parameters = [number_of_spots, pixel_size, pixel_size*image_size[0], image_size[0], wavelength, NA, objectiveFocalLength, [n_sample, n_objective, n_objective], nDiscretizationBFP, background, n_photons, x_pos_nm, y_pos_nm, theta, phi]

    # Return results instead of printing
    return dipole_psf_noisy, parameters


def run_mortensen_fit(image, image_path, parameters_path): # change this to accept image + parameters files

    psf_image_output = tifffile.imread(image_path)
#    psf_image_output = np.load(image_path)
    psf_image_simulated = image

    psf_image = psf_image_output

    parameters = parse_parameters_file(parameters_path)
    theta = parameters['angleInclination_array'][0]
    phi =  parameters['angleAzimuth_array'][0]
    x_pos_nm = parameters['positionX_nm_array'][0]
    y_pos_nm = parameters['positionY_nm_array'][0]
    n_photons = parameters['par.nPhotons']

    image_size = (parameters['image_size_px'], parameters['image_size_px'])
    pixel_size = parameters['pixel_size_nm']
    wavelength = parameters['wavelength']
    n_objective = parameters['par.refractiveIndices'][1]
    n_sample = parameters['par.refractiveIndices'][0]
    magnification = 215
    NA = parameters['par.objectiveNA']

    norm_file = "/home/tfq96423/dipolenorm.npy"

    # factor by which to oversample
    subpixel_factor = 1

#    print(image_size)
#    print(pixel_size)
#    print(wavelength)
#    print(n_objective)
#    print(n_sample)
#    print(magnification)
#    print(NA)
#    print(n_photons)



    print('-------------------')
    print(f'ϕ = {round(phi*180/np.pi)}°/360°')

    start_time = time.time()

    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor)    

    # Run Mortensen fit
    results = psf_generator.mortensen_fit(psf_image)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'    {elapsed_time:.2f} seconds')

    
    return results, [phi, theta, x_pos_nm, y_pos_nm, n_photons]



# Prevent script from running when imported
if __name__ == "__main__":

    runs = range(1)
    thetas = np.arange(0, np.pi/2 + 0.001, 10*np.pi/180)#np.arange(0, np.pi/2 + 0.001, 22.5*np.pi/180)
    phis = [0]#np.arange(0, 2*np.pi - 0.001, 45*np.pi/180)

    for run in runs:
      for theta in thetas:
        for phi in phis:

          theta_deg = theta*180/np.pi
          phi_deg = phi*180/np.pi

          # Define output paths
          image_output_path = f'/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/animation_highres/high_res_images/mortensen_theta{ceil(theta_deg):03}.tif'
          params_output_path = f"/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/animation_highres/high_res_images/params_theta{ceil(theta_deg):03}_phi{ceil(phi_deg):03}_run{int(run+1)}.m"

          psf_image, parameters = run_simulator(phi, theta)
          save_as_tif(psf_image, image_output_path)
          # save the array itself
#          np.save(image_output_path, psf_image)    
#          np.savetxt(image_output_path, psf_image, delimiter=",")

          # Write data to file
          with open(params_output_path, 'w') as file:
            file.write(f"% ground truth for sim_inc{round(theta_deg)}_az{round(phi_deg)}_run{int(run+1)}.tif\n")
            file.write("% settings\n")
            file.write(f"number_of_spots = {parameters[0]}\n")
            file.write(f"pixel_size_nm = {parameters[1]}\n")
            file.write(f"image_size_nm = {parameters[2]}\n")
            file.write(f"image_size_px = {parameters[3]}\n")
            file.write(f"wavelength = {parameters[4]}\n")
            file.write(f"par.objectiveNA = {parameters[5]}\n")
            file.write(f"objectiveFocalLength = {parameters[6]}\n")
            file.write(f"par.refractiveIndices = {parameters[7]}\n")
            file.write(f"par.nDiscretizationBFP = {parameters[8]}\n")
            file.write(f"par.backgroundNoise = {parameters[9]}\n")
            file.write(f"par.nPhotons = {parameters[10]}\n")
            file.write(f"positionX_nm_array = [{parameters[11]}, ]\n")
            file.write(f"positionY_nm_array = [{parameters[12]}, ]\n")
            file.write(f"angleInclination_array = [{parameters[13]}, ]\n")
            file.write(f"angleAzimuth_array = [{parameters[14]}, ]\n")
          file.close()

