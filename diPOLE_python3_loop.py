'''
This script uses Mortensen's dipole fitting code diPOLE.py from 'Supplementary Software' at the bottom of https://www.nature.com/articles/ncomms9621
I've translated it into python3 to work, and modified Estimate() so that it returns the position, orientation, and covariance matrix.
That code is for a single spot on a single frame, so here we just loop over every blob in a frame, and then over every frame.

One complication is that to do this we need to be able to identify spots in a frame. For this, we rely on thunderSTORM.
This script runs thunderSTORM on a directory of frames. We take those localisations, extract an NxN pixel region centred on them, and then apply Mortensen to that patch.
We then replace the localisations in the thunderSTORM results with these new adjusted localisations.
The output is a typical thunderSTORM results table. I tried to automatically turn this into an image with imageJ, but couldn't.

One of the main issues currently is that I've struggled to integrate ImageJ here.
It runs thunderSTORM fine, but I would prefer it to run headless if possible.
I also tried to get it to run at the end to generate the visualisation of the results, but have been hitting a wall with that.
pyimagej seems difficult to work with, but it might just be me.

Also, I don't have the right parameters for the experimental setup,
or for the initial thunderSTORM run. So correct that.
And make it read them in from somewhere.
Need to read in the pixel scale from the image metadata or something. Is that available?

Running on /mnt/cryosil_ro/AttoDRY800/Developments_Autumn_2024/2024-10-09_ASIL240923C05/StormData/StormData_1/Results
and using the parameters from protocol.txt

Oh and it's super slow. Like 15 seconds per spot / 5 minutes per frame. So this would never work for our usual 10,000 frame stack

--------------------

First: run thunderstorm on all frames. This gets the initial localisations.
Then loop over every frame i:
    Run mortensen_single_frame()
      This will consider all x,y in the thunderstorm results which have frame=i
      It will run extract_patches() on frame=i
      It will loop over all patches in that frame, run Mortensen on it
      Append the resulting (framenumber, xs, ys) to the overall array
Add the xs and ys to the existing thunderstorm results table
(because Mortensen is done relative to a small patch centred on the thunderstorm localisations, so just need to add it on)
'''

import diPOLE_python3
from pylab import *
import numpy as np
import cv2
import pandas as pd
import subprocess
import os
import time

# from config_params import *
from thunderstorm_run_macro import run_thunderstorm, reconstruct
# from thunderstorm_reconstruct_macro import reconstruct

def blob_detect_all_frames(frames_dir, results_path, pixel_width):

    fs, xs, ys = run_thunderstorm(frames_dir, results_path)

    # rescale from nm to px
    xs_image_coords = xs/pixel_width
    ys_image_coords = ys/pixel_width

    centroids_image_coords = [(int(f), int(x), int(y)) for f, x, y in zip(fs, xs_image_coords, ys_image_coords)]

    return centroids_image_coords


def extract_patches(image, current_frame_number, centroids_image_coords, patch_width):
    """
    Extract patches around given centroids from the image.

    Args:
        image (ndarray): The input image from which to extract the patches.
        centroids (list of tuple): A list of (x, y) coordinates of the centroids.
        patch_size (int): The size of the patch to extract (default is 12).

    Returns:
        list of ndarray: A list of extracted patches, or None for out of bounds.
    """
    patches = []

    for centroid in centroids_image_coords:
        f, x, y = centroid

        # only operate on the current frame number
        if f == current_frame_number:

            # Define the coordinates for the patch
            x_start = x - patch_width // 2
            x_end = x + patch_width // 2
            y_start = y - patch_width // 2
            y_end = y + patch_width // 2

            # # if you decide to ignore out-of-bounds patches:
            # # Check for out-of-bounds
            # if x_start >= 0 and x_end <= image.shape[1] and y_start >= 0 and y_end <= image.shape[0]:
            #     # Extract the patch
            #     patch = image[y_start:y_end, x_start:x_end]
            #     patches.append(patch)

            # if you decide to pad out-of-bounds patches:
            # Determine the padding needed for each side
            pad_left = max(0, -x_start)  # How many pixels need padding on the left
            pad_right = max(0, x_end - image.shape[1])  # How many pixels need padding on the right
            pad_top = max(0, -y_start)  # How many pixels need padding on the top
            pad_bottom = max(0, y_end - image.shape[0])  # How many pixels need padding on the bottom

            # Apply padding to the image
            padded_image = np.pad(image,
                                  ((pad_top, pad_bottom), (pad_left, pad_right)),
                                  mode='constant', constant_values=0)

            # Recalculate the new patch coordinates based on the padded image
            x_start_padded = x_start + pad_left
            x_end_padded = x_end + pad_left
            y_start_padded = y_start + pad_top
            y_end_padded = y_end + pad_top

            # Extract the patch from the padded image
            patch = padded_image[y_start_padded:y_end_padded, x_start_padded:x_end_padded]
            patches.append(patch)

    return patches


def mortensen_single_frame(image_path,
                           current_frame_number,
                           centroids_image_coords,
                           patch_width,
                           peak_emission_wavelength,
                           pixel_width,
                           magnification,
                           numerical_aperture,
                           ref_ind_immersion,
                           ref_ind_imaging,
                           ref_ind_buffer,
                           initvals,
                           initpix,
                           deltapix,
                           Sfloor,
                           inverse_gain,
                           sigma_noise):

    # Load image
    image = cv2.imread(image_path, 0)

    # Subtract the background roughly
    image = np.clip(image - np.mean(image), 0, 255).astype(np.uint8)

    # get 12x12 patches around each centroid
    blob_patches = extract_patches(image, current_frame_number, centroids_image_coords, patch_width)

    # Create instance of MLEwT !!! this was inside the loop earlier !!!
    track = diPOLE_python3.MLEwT(peak_emission_wavelength,
                                 pixel_width,
                                 magnification,
                                 numerical_aperture,
                                 ref_ind_immersion,
                                 ref_ind_imaging,
                                 ref_ind_buffer,
                                 initvals,
                                 initpix,
                                 deltapix,
                                 Sfloor,
                                 inverse_gain,
                                 sigma_noise)

    # Loop over all blobs
    x_list, y_list, phi_list, theta_list, covariance_list = [], [], [], [], []

    # for blob in blob_patches:
    for i, blob in enumerate(blob_patches, 1):

        start_blob = time.time()
        print(f"Analysing blob {i}/{len(blob_patches)}")

        # Perform estimation
        x_est, y_est, theta_est, phi_est, cov_mat = track.Estimate(blob)
        x_list.append(x_est)
        y_list.append(y_est)
        theta_list.append(theta_est)
        phi_list.append(phi_est)
        covariance_list.append(cov_mat)

        end_blob = time.time()
        elapsed_time_blob = end_blob - start_blob
        print(f"Time: {elapsed_time_blob:.4f} seconds on this blob")

    return x_list, y_list, theta_list, phi_list, covariance_list


# --------------------
# Mortensen run params
# --------------------
# Experimental parameters
peak_emission_wavelength = 580.0 # Peak emission wavelength
ref_ind_immersion = 1.52 # RI of immersion oil/optics
ref_ind_imaging = 1.0 # Refractive index of imaging medium (typically that of air, i.e. np = 1.0)
ref_ind_buffer = 1.33 # RI of buffer
numerical_aperture = 1.49 # Numerical aperture of objective
magnification = 250.0 # Magnification of composite microscope
pixel_width = 51.2 # Pixel width (nm per px)

# PSF parameters
photon_number = 10000.0 # Photon number?
background_level = 5.0 # Background level?
mu = 0.1 # Probe location?
nu = 0.1 # ...
phi = 2 * np.pi / 3.0 # inclination
theta = 0.5 # azimuthal angle
deltaz = -30.0 # Distance from design focal plane

# EMCCD parameters
inverse_gain = 1./100.
sigma_noise = 2 #12.6
Sfloor = 300.0
gain = 1.0 / inverse_gain

patch_width = 12 # size of NxN pixel patch around blob centroid to consider
# --------------------

# Initial thunderstorm run to get blob location

# Load data
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/ASIL240923C05_50.ome.tif-frames/'
results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/thunderstorm_results.csv'

# find centroids using gaussian fitting thunderstorm
centroids_image_coords = blob_detect_all_frames(frames_dir, results_path, pixel_width)

# get frames
frame_paths = [os.path.join(frames_dir, f) for f in os.listdir(frames_dir) if f.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'))]

# Initial guess params
initvals = array([mu, nu, background_level, photon_number, phi, theta, deltaz]) # initial PSF values
deltapix = patch_width / 2 # centre of patch around blob
initpix = (deltapix, deltapix) # centre of patch around blob

# Mortensen run on each blob in each frame
x_ests, y_ests, theta_ests, phi_ests, covariance_ests = [], [], [], [], []
for i, frame_path in enumerate(frame_paths, 1):

    print("----------")
    print(f"FRAME {i}")
    start_frame = time.time()

    single_frame_results = list(mortensen_single_frame(frame_path,
                               i,
                               centroids_image_coords,
                               patch_width,
                               peak_emission_wavelength,
                               pixel_width,
                               magnification,
                               numerical_aperture,
                               ref_ind_immersion,
                               ref_ind_imaging,
                               ref_ind_buffer,
                               initvals,
                               initpix,
                               deltapix,
                               Sfloor,
                               inverse_gain,
                               sigma_noise))

    x_ests.append(single_frame_results[0])
    y_ests.append(single_frame_results[1])
    theta_ests.append(single_frame_results[2])
    phi_ests.append(single_frame_results[3])
    covariance_ests.append(single_frame_results[4])

    end_frame = time.time()
    elapsed_time_frame = end_frame - start_frame
    elapsed_time_frame = elapsed_time_frame/60
    print(f"Time: {elapsed_time_frame:.4f} minutes on this frame")

# make sure list is flat, because it needs to be for results table
x_ests = [item for sublist in x_ests for item in sublist]
y_ests = [item for sublist in y_ests for item in sublist]
theta_ests = [item for sublist in theta_ests for item in sublist]
phi_ests = [item for sublist in phi_ests for item in sublist]
covariance_ests = [item for sublist in covariance_ests for item in sublist]

print(len(x_ests))
print(len(y_ests))
print(len(theta_ests))
print(len(phi_ests))
print(len(covariance_ests))

# generate a thunderstorm-style results table with new x,y localisations
mortensen_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/mortensen_results.csv'
# output_img_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/reconstruction.png'

df = pd.read_csv(results_path)
if len(x_ests) != len(df) or len(y_ests) != len(df):
    raise ValueError("The length of the new x and y arrays must match the number of rows in the CSV file.")
df['x [nm]'] += x_ests
df['y [nm]'] += y_ests
df.to_csv(mortensen_results_path, index=False)

# --------------------
# generating image from results table
# [not working - go and do it in fiji manually for now]
# --------------------
# # this doesn't work due to some conflict or something
# reconstruct(mortensen_results_path, output_img_path)
# # try getting round it by running the macro via terminal
# # !!! this has input/output paths hard-coded !!!
# command = "/home/tfq96423/fiji-linux64/Fiji.app/ImageJ-linux64 -macro /home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/reconstruct.ijm"
# subprocess.run(command, shell=True, check=True)
#
# # show in napari
# command2 = f"napari {output_img_path}"
# subprocess.run(command2, shell=True, check=True)
