import numpy as np
import os
import sys
import time
import re
import datetime
import matplotlib.pyplot as plt
from MLEwT_fixed import dipdistr, MLEwT
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

class DipolePSFGenerator:
    def __init__(self, image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, verbose=False):
        self.image_size = image_size
        self.pixel_size = pixel_size
        self.wavelength = wavelength
        self.n_sample = n_sample
        self.n_objective = n_objective
        self.magnification = magnification
        self.NA = NA
        self.norm_file = norm_file
        self.verbose = verbose
        
        # Load normalization file
        self.norm_data = np.load(norm_file)
        
        # Initialize dipole distribution model
        self.DD = dipdistr(wavelength, n_objective, n_sample, magnification, NA, norm_file)
    
    def __call__(self, phi, theta, x_pos, y_pos, n_photons):
        
        # Create position vector
        # x_pos and y_pos are in nm
        posvec = np.arange(-(self.image_size[0]-1)/2, self.image_size[0]/2) * self.pixel_size
        dipole_psf = np.zeros(self.image_size)
        
        # Generate PSF
        for i in range(self.image_size[0]):
            for j in range(self.image_size[1]):
                dipole_psf[j, i] = self.DD.PSF_approx(posvec[i] - x_pos, 
                                                      posvec[j] - y_pos,
                                                      phi, theta, 
                                                      )
        
        dipole_psf = dipole_psf / dipole_psf.sum() * n_photons
        
        return dipole_psf

    def mortensen_fit(self, dipole_psf):
        
        # Define parameters
        #deltapix = 9
        
        # Generate the initial mux, muy in pixels around the center of the image (9, 9)
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
        #print(f"initial values for the center pixel (ypixel,xpixel): {initpix}")
        
        # initvals = the initial values to start the estimate       
        track = MLEwT(initvals, self)   #initpix, deltapix, self)
        result = track.Estimate(dipole_psf)

        return result

def run_mortensen_fit(image_path, parameters_path): # change this to accept image + parameters files

    psf_image = tifffile.imread(image_path)

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


    print(image_size)
    print(pixel_size)
    print(wavelength)
    print(n_objective)
    print(n_sample)
    print(magnification)
    print(NA)
    print(n_photons)



    print('-------------------')
    print(f'ϕ = {round(phi*180/np.pi)}°/360°')

    start_time = time.time()

    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file)    

    # Run Mortensen fit
    results = psf_generator.mortensen_fit(psf_image)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'    {elapsed_time:.2f} seconds')

    
    return results, [phi, theta, x_pos_nm, y_pos_nm, n_photons]

# Prevent script from running when imported
if __name__ == "__main__":

    # find each tif and param file in sims directory
    # loop over each one
    # run fitting process on each one    

    #results_path = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_mortensen/fitting_results_mortensen.py"
    results_path = "./fitting_results_mortensen.py"

    # find images and parameter files

    valid_extensions_image = {'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'}
    valid_extensions_settings = {'.m'}

    frame_paths = []
    settings_paths = []

    frames_dir = './sims'  # Replace with your actual directory

    # Iterate over files in the directory
    for filename in os.listdir(frames_dir):
        # Get the file extension
        _, ext = os.path.splitext(filename)
        ext = ext.lower()
    
        # If the file is an image, add its path to frame_paths
        if ext in valid_extensions_image:
            frame_paths.append(os.path.join(frames_dir, filename))
    
        # If the file is a settings file, add its path to settings_paths
        elif ext in valid_extensions_settings:
            settings_paths.append(os.path.join(frames_dir, filename))



    # run the fits on every frame
    for i in range(len(frame_paths)):
      results, ground_truth = run_mortensen_fit(frame_paths[i], settings_paths[i])
      save_results_to_py(results, ground_truth, results_path)



