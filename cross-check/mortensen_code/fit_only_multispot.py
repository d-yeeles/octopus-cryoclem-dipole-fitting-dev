import numpy as np
from MLEwT_fixed import dipdistr, MLEwT
from math import ceil, floor
import os
import sys
import time
import re
import datetime
import matplotlib.pyplot as plt
import tifffile
from save_results import save_results_to_py

def normalize_image(image):
    """
    Normalize image exactly like the MATLAB code:
    (image - min) / (max - min)
    """
    # Convert to float64 for calculations
    img_float = image.astype(np.float64)

    # Find min and max values
    min_val = np.min(img_float)
    max_val = np.max(img_float)

    # Avoid division by zero
    if max_val == min_val:
        return np.zeros_like(img_float)

    # Normalize to 0-1 range exactly like the MATLAB code
    normalized = (img_float - min_val) / (max_val - min_val)

    # Convert to 8-bit (0-255)
    img_8bit = (normalized * 255).astype(np.uint8)

    return img_8bit

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





    def mortensen_fit(self, initpix, deltapix, dipole_psf):
  
        # Generate the initial mux, muy in pixels around the centre of the image (9, 9)
#        mux_pix = np.random.uniform(self.image_size[1] / 2 - 0.5, self.image_size[1]/ 2 + 0.5)
#        muy_pix = np.random.uniform((self.image_size[0] - 1)/ 2 - 0.5, (self.image_size[0] - 1) / 2 + 0.5)
        #print(f"the initial vals of mux muy: {mux_pix, muy_pix}")
        mux_pix, muy_pix = initpix

        # Convert to coordinates in nm
#        mux_nm = np.random.uniform(0 - self.pixel_size/2, 0 + self.pixel_size/2)    
#        muy_nm = np.random.uniform(0 - self.pixel_size/2, 0 + self.pixel_size/2)    
        mux_nm = (mux_pix - (self.image_size[1] - 1) / 2) * self.pixel_size
        muy_nm = (muy_pix - (self.image_size[0] - 1) / 2) * self.pixel_size

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
        track = MLEwT(initvals, initpix, deltapix, self)
        result = track.Estimate(dipole_psf)

        return result

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

    parameters = [number_of_spots, pixel_size, pixel_size*image_size[0], image_size[0], wavelength, NA, objectiveFocalLength, [n_sample, n_objective, n_objective], nDiscretizationBFP, background, n_photons, x_pos_nm_global, y_pos_nm_global, theta, phi]

    # Return results instead of printing
    return dipole_psf_noisy, parameters


def run_mortensen_fit(initpix, deltapix, full_image, frame_parameters, dipole_index):
    """
    Run the Mortensen fit on a patch extracted from the full image.
    
    Parameters:
    -----------
    initpix : tuple (y, x)
        Initial pixel coordinates in the full image
    deltapix : int
        Half-size of the patch to extract
    full_image : ndarray
        The full image containing all dipoles
    frame_parameters : dict
        Parameters loaded from the .m file
    dipole_index : int
        Index of the dipole being analyzed
    """
    # Read ground truth from params file for comparison
    theta = frame_parameters['angleInclination_array'][dipole_index]
    phi = frame_parameters['angleAzimuth_array'][dipole_index]
    x_pos_nm = frame_parameters['positionX_nm_array'][dipole_index]
    y_pos_nm = frame_parameters['positionY_nm_array'][dipole_index]
    n_photons = frame_parameters['par.nPhotons']
    
    # Get necessary parameters
    image_size = (frame_parameters['image_size_px'], frame_parameters['image_size_px'])
    pixel_size = frame_parameters['pixel_size_nm']
    wavelength = frame_parameters['wavelength']
    n_objective = frame_parameters['par.refractiveIndices'][1]
    n_sample = frame_parameters['par.refractiveIndices'][0]
    magnification = 215
    NA = frame_parameters['par.objectiveNA']
    norm_file = "/home/tfq96423/dipolenorm.npy"
    subpixel_factor = 1
    
    # Calculate image centre in pixels
    full_centre_y = full_image.shape[0] / 2
    full_centre_x = full_image.shape[1] / 2
    
    # Extract the patch with bounds checking
    ypix, xpix = initpix

#    print(initpix)
#    print(deltapix)

    y_start = max(0, ypix - deltapix)
    y_end = min(full_image.shape[0], ypix + deltapix+1)
    x_start = max(0, xpix - deltapix)
    x_end = min(full_image.shape[1], xpix + deltapix+1)
    
    # Extract the patch image
    patch_image = full_image[y_start:y_end, x_start:x_end]
    patch_height, patch_width = patch_image.shape

    # Need to compare the patch edges and centre to the whole image - they should be identical here
#    print(full_image.shape)
#    print(patch_image.shape)
    
    # Define patch centre in PATCH coordinates
    patch_centre_y = patch_height / 2  
    patch_centre_x = patch_width / 2
    
    # Create PSF generator with PATCH dimensions
    patch_size = (patch_height, patch_width)
    psf_generator = DipolePSFGenerator(
        patch_size, pixel_size, wavelength, n_objective, n_sample, 
        magnification, NA, norm_file, subpixel_factor
    )

    # Initialize fitting parameters
    # Random initialization for orientation parameters
    init_new1 = np.random.uniform(-1, 1)
    init_new2 = np.random.uniform(-1, 1)
    init_new3 = np.random.uniform(0, 1)
    
    # Initialize position at patch centre
    mux_nm = 0.0
    muy_nm = 0.0
    
    # Estimate initial photon count
    init_photons = np.sum(patch_image)
    
    # Create initial values array for MLEwT
    initvals = np.array([init_new1, init_new2, init_new3, mux_nm, muy_nm, init_photons])
    
    # Center pixel in the patch coordinates
    patch_initpix = (floor(patch_centre_y), floor(patch_centre_x))
   
#    # Ensure patch_initpix is within patch bounds
#    patch_initpix = (
#        min(max(0, patch_initpix[0]), patch_height - 1),
#        min(max(0, patch_initpix[1]), patch_width - 1)
#    )

    # Run the MLEwT fit (which returns local coordinates in nm relative to patch centre)
#    track = MLEwT(initvals, patch_initpix, min(deltapix, patch_height//2, patch_width//2), psf_generator)
    track = MLEwT(initvals, patch_initpix, deltapix, psf_generator)
    local_results = track.Estimate(patch_image)
    
    # Get fitted positions in nm relative to patch centre
    fitted_x_nm_local = local_results[2]
    fitted_y_nm_local = local_results[3]
    
    # ---------- COORDINATE TRANSFORMATION SECTION ----------
    
    # 1. Calculate position of patch centre in global image coordinates (pixels)
    patch_centre_global_y = y_start + patch_centre_y
    patch_centre_global_x = x_start + patch_centre_x
    
    # 2. Calculate position of patch centre in global coordinates (nm)

#    print('patch_centre_global px', patch_centre_global_y, patch_centre_global_x)
#    print('full_centre px', full_centre_y, full_centre_x)

    patch_centre_global_y_nm = (patch_centre_global_y - full_centre_y) * pixel_size
    patch_centre_global_x_nm = (patch_centre_global_x - full_centre_x) * pixel_size
    
    # 3. Add local offsets (from fitting) to get final global positions
    fitted_y_nm_global = patch_centre_global_y_nm + fitted_y_nm_local
    fitted_x_nm_global = patch_centre_global_x_nm + fitted_x_nm_local
    
    # Create final results with corrected positions
    global_results = local_results.copy()
    global_results[2] = fitted_x_nm_global
    global_results[3] = fitted_y_nm_global
    
#    # Print debugging info
#    print(f"POSITION TRACKING:")
#    print(f"  Patch extraction: y=[{y_start}:{y_end}], x=[{x_start}:{x_end}]")
#    print(f"  Patch dimensions: {patch_image.shape}")
#    print(f"  Initial pixel in full image: {initpix}")
#    print(f"  Initial pixel in patch: {patch_initpix}")
#    print(f"  Patch centre (global px): y={patch_centre_global_y}, x={patch_centre_global_x}")
#    print(f"  Patch centre (global nm): y={patch_centre_global_y_nm}, x={patch_centre_global_x_nm}")
#    print(f"  Ground truth global (nm): x={x_pos_nm}, y={y_pos_nm}")
#    print(f"  Fitted local position (nm): x={fitted_x_nm_local}, y={fitted_y_nm_local}")
#    print(f"  Transformed to global (nm): x={fitted_x_nm_global}, y={fitted_y_nm_global}")
#    print(f"  Full image dimensions (px): {full_image.shape}")
#    print(f"  Pixel size (nm): {pixel_size}")
   

#    # Inside run_mortensen_fit
#    plt.figure(figsize=(10, 5))
#    plt.subplot(1, 2, 1)
#    plt.imshow(patch_image, cmap='viridis')
#    plt.title(f"Patch at ({initpix[0]}, {initpix[1]})")
#    plt.colorbar()
#    
#    # Also visualize the fitted PSF
#    fitted_psf = psf_generator(
#        global_results[0],  # phi
#        global_results[1],  # theta
#        0,  # x_pos (patch centre)
#        0,  # y_pos (patch centre)
#        global_results[4]   # n_photons
#    )
#    plt.subplot(1, 2, 2)
#    plt.imshow(fitted_psf, cmap='viridis')
#    plt.title("Fitted PSF")
#    plt.colorbar()
#    plt.savefig(f"patch_debug_{dipole_index}.png")
#    plt.close()


 
    return global_results, [phi, theta, x_pos_nm, y_pos_nm, n_photons]


# Prevent script from running when imported
if __name__ == "__main__":

    # For simulated tests:
    # --------------------
    # This does the following:
    #  - Take a dir full of frames and their associated params
    #  - For each frame in the dir:
    #     - Read in frame, params
    #     - For each position in positionX_nm (or whatever)
    #        - Set that as initpix (+/- some simulated error)
    #        - Run the fitting
    #        - Output results

    # For real images:
    # --------------------
    # Rewrite this so that it does the following:
    #  - Take a tif stack
    #  - For each frame in the stack:
    #     - Read in frame, run thunderstorm to get localisations
    #     - For each thunderstorm localisation:
    #        - Set that as initpix
    #        - Run the fitting
    #        - Output results

    # Output path
    results_path = '../results_1spot_loads_normalised/results_1spot_mortensen_normalised.py'

    # find images and parameter files
    frames_dir = '../sims_1spot_loads'

    valid_extensions_image = {'.tif'}
    valid_extensions_settings = {'.m'}

    frame_paths = []
    settings_paths = []

    # Iterate over files in the directory
    for filename in sorted(os.listdir(frames_dir)):
        # Get the file extension
        _, ext = os.path.splitext(filename)
        ext = ext.lower()
    
        # If the file is an image, add its path to frame_paths
        if ext in valid_extensions_image:
            frame_paths.append(os.path.join(frames_dir, filename))
    
        # If the file is a settings file, add its path to settings_paths
        elif ext in valid_extensions_settings:
            settings_paths.append(os.path.join(frames_dir, filename))

    # Loop over each frame in dir
    for frame_index in range(len(frame_paths)):

        print('-------------------')
        print(f'FRAME {frame_index+1}/{len(frame_paths)}')

        # Read in image
        psf_image = tifffile.imread(frame_paths[frame_index])
        #psf_image = np.loadtxt(frame_paths[frame_index], delimiter=",")

        # Read in parameters
        frame_parameters = parse_parameters_file(settings_paths[frame_index])
#        print(frame_parameters['image_size_nm'])
        for_display_xest_list = []
        for_display_yest_list = []
        for_display_xtru_list = []
        for_display_ytru_list = []

        # Loop over each dipole in frame
        for dipole_index in range(len(frame_parameters['positionX_nm_array'])):

            # Use ground truth position + artificial error as initpix

            # Calculate global position in pixels (relative to the centre of the grid)
            grid_size_px = frame_parameters['image_size_px']
            grid_centre_x = psf_image.shape[1] / 2
            grid_centre_y = psf_image.shape[0] / 2

            # Convert from nm (relative to grid centre) to pixels (relative to top-left)
            initpix_x = (frame_parameters['positionX_nm_array'][dipole_index] / frame_parameters['pixel_size_nm']) + grid_centre_x
            initpix_y = (frame_parameters['positionY_nm_array'][dipole_index] / frame_parameters['pixel_size_nm']) + grid_centre_y

            initpix = (int(initpix_y), int(initpix_x))  # Note: (y,x) ordering for numpy indexing






            # Run the fit
            print(f"      {dipole_index+1}/{len(frame_parameters['positionX_nm_array'])}")
            start_time = time.time()
            deltapix = 9
            results, ground_truth = run_mortensen_fit(initpix, deltapix, psf_image, frame_parameters, dipole_index)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f'    {elapsed_time:.2f} seconds')

            # Output results
            save_results_to_py(results, ground_truth, results_path)

#            for_display_xest_list.append(results[2])
#            for_display_yest_list.append(results[3])
#            for_display_xtru_list.append(ground_truth[2])
#            for_display_ytru_list.append(ground_truth[3])
#
#
#
#        # Diagnostic plot
#        display_psf_image = normalize_image(psf_image)
#        plt.figure(figsize=(24, 24))
#        plt.imshow(display_psf_image, cmap='gray', vmin=0, vmax=255, extent=[-frame_parameters['image_size_nm']/2, frame_parameters['image_size_nm']/2, -frame_parameters['image_size_nm']/2, frame_parameters['image_size_nm']/2])
#        plt.scatter(for_display_xtru_list, for_display_ytru_list, c='green')
#        plt.scatter(for_display_xest_list, for_display_yest_list, c='red')
#        plt.savefig(f"../results_1ispot/diagnostic_frame{frame_index+1}.png")
#        plt.close()




