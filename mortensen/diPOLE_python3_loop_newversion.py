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

# # from diPOLE_python3_loop import centroids_image_coords
# from config_params import *
# from thunderstorm_run_macro import run_thunderstorm, reconstruct
# # from thunderstorm_reconstruct_macro import reconstruct

def nm_to_px(coord_nm, pixel_size_nm, image_size_px, x_or_y_flag):
    if x_or_y_flag == 'x':
        coord_px = (coord_nm / pixel_size_nm + (image_size_px + 1) / 2) - 1
    elif x_or_y_flag == 'y':
        coord_px = image_size_px - (coord_nm / pixel_size_nm + image_size_px / 2) - 1
    return coord_px

# def blob_detect_all_frames(frames_dir, results_path, pixel_size_nm):
#
#     run_thunderstorm(frames_dir, results_path)
#
#     # Import the positions in each frame
#     df = pd.read_csv(results_path)
#     df_selected_columns = df[['frame','x [nm]', 'y [nm]']]
#     f_array = df['frame'].to_numpy()
#     x_array = df['x [nm]'].to_numpy()
#     y_array = df['y [nm]'].to_numpy()
#
#     # rescale from nm to px
#     x_array_image_coords = x_array/pixel_size_nm
#     y_array_image_coords = y_array/pixel_size_nm
#
#     centroids_image_coords = [(int(f), int(x), int(y)) for f, x, y in zip(f_array, x_array_image_coords, y_array_image_coords)]
#
#     return centroids_image_coords


# def replace_outside_patch_with_zeros(image, current_frame_number, centroids_nm, patch_width_px):
#     """
#     Replace pixels outside the defined patch area with 0s and return the masked images as a list.
#
#     Args:
#         image (ndarray): The input image.
#         centroids_nm (list of tuple): A list of (f, x, y) where f is the frame number, and (x, y) are the coordinates of the centroids.
#         patch_width_px (int): The width of the patch.
#         pixel_size_nm (float): The pixel size in nm.
#
#     Returns:
#         list of ndarray: A list of masked images where pixels outside the patch area are replaced with 0s.
#     """
#     masked_images = []  # List to store the masked images
#     image_with_rectangles = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)
#
#     # Iterate over the centroids
#     for centroid in centroids_nm:
#         f, x, y = centroid
#
#         # Only operate on the current frame number
#         if f == current_frame_number:
#
#             # Convert the centroids from nm to px
#             x_px = int(nm_to_px(x, pixel_size_nm, image.shape[1], 'x'))
#             y_px = int(nm_to_px(y, pixel_size_nm, image.shape[1], 'y'))
#
#             # Define the coordinates for the patch
#             x_start = x_px - patch_width_px // 2
#             x_end = x_px + patch_width_px // 2
#             y_start = y_px - patch_width_px // 2
#             y_end = y_px + patch_width_px // 2
#
#             # Create a mask of the region inside the patch
#             mask = np.zeros_like(image, dtype=np.uint8)
#
#             # Make sure the mask covers the patch area
#             mask[y_start:y_end, x_start:x_end] = 1
#
#             # Apply the mask to the image (set the pixels outside the patch to 0)
#             image_masked = image * mask
#
#             # Add the masked image to the list
#             masked_images.append(image_masked)
#
#             # Optional: Draw a rectangle around the patch on the image
#             cv2.rectangle(image_with_rectangles, (x_start, y_start), (x_end, y_end), (255, 0, 0), 1)
#
#     # Optional: Display the image with the rectangles
#     x_coords_nm = [centroid[1] for centroid in centroids_nm]
#     y_coords_nm = [centroid[2] for centroid in centroids_nm]
#     centroids_x_px = np.array([nm_to_px(x, pixel_size_nm, image.shape[1], 'x') for x in x_coords_nm])
#     centroids_y_px = np.array([nm_to_px(y, pixel_size_nm, image.shape[1], 'x') for y in y_coords_nm])
#
#     fig, ax = plt.subplots(figsize=(5, 5))
#     ax.imshow(image_with_rectangles, cmap='gray')
#     ax.scatter(centroids_x_px, centroids_y_px, color='red', s=10)  # s is the size of the dot
#     ax.axis('off')
#     plt.show()
#
#     return masked_images

# def extract_patches(image, current_frame_number, centroids_nm, patch_width_px):
#     """
#     Extract patches around given centroids from the image.
#
#     Args:
#         image (ndarray): The input image from which to extract the patches.
#         centroids (list of tuple): A list of (x, y) coordinates of the centroids.
#         patch_size (int): The size of the patch to extract (default is 12).
#
#     Returns:
#         list of ndarray: A list of extracted patches, or None for out of bounds.
#     """
#     patches = []
#
#     image_with_rectangles = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)
#
#     for centroid in centroids_nm:
#         f, x, y = centroid
#
#         # only operate on the current frame number
#         if f == current_frame_number:
#
#             x_px = int(nm_to_px(x, pixel_size_nm, image.shape[1], 'x'))
#             y_px = int(nm_to_px(y, pixel_size_nm, image.shape[1], 'y'))
#
#             # Define the coordinates for the patch
#             x_start = x_px - patch_width_px // 2
#             x_end = x_px + patch_width_px // 2
#             y_start = y_px - patch_width_px // 2
#             y_end = y_px + patch_width_px // 2
#
#             # # if you decide to ignore out-of-bounds patches:
#             # # Check for out-of-bounds
#             # if x_start >= 0 and x_end <= image.shape[1] and y_start >= 0 and y_end <= image.shape[0]:
#             #     # Extract the patch
#             #     patch = image[y_start:y_end, x_start:x_end]
#             #     patches.append(patch)
#
#             # if you decide to pad out-of-bounds patches:
#             # Determine the padding needed for each side
#             pad_left = max(0, -x_start)  # How many pixels need padding on the left
#             pad_right = max(0, x_end - image.shape[1])  # How many pixels need padding on the right
#             pad_top = max(0, -y_start)  # How many pixels need padding on the top
#             pad_bottom = max(0, y_end - image.shape[0])  # How many pixels need padding on the bottom
#
#             # Apply padding to the image
#             padded_image = np.pad(image,
#                                   ((pad_top, pad_bottom), (pad_left, pad_right)),
#                                   mode='constant', constant_values=0)
#
#             # Recalculate the new patch coordinates based on the padded image
#             x_start_padded = x_start + pad_left
#             x_end_padded = x_end + pad_left
#             y_start_padded = y_start + pad_top
#             y_end_padded = y_end + pad_top
#
#             # Extract the patch from the padded image
#             patch = padded_image[y_start_padded:y_end_padded, x_start_padded:x_end_padded]
#             patches.append(patch)
#
#             # uncomment to display patches
#             cv2.rectangle(image_with_rectangles, (x_start_padded, y_start_padded), (x_end_padded, y_end_padded), (255, 0, 0), 1)
#     x_coords_nm = [centroid[1] for centroid in centroids_nm]
#     y_coords_nm = [centroid[2] for centroid in centroids_nm]
#     centroids_x_px = np.array([nm_to_px(x, pixel_size_nm, image.shape[1], 'x') for x in x_coords_nm])
#     centroids_y_px = np.array([nm_to_px(y, pixel_size_nm, image.shape[1], 'x') for y in y_coords_nm]) # using 'x' because scatter uses bottom-left origin
#     fig, ax = plt.subplots(figsize=(5, 5))
#     ax.imshow(image_with_rectangles, cmap='gray')
#     ax.scatter(centroids_x_px, centroids_y_px, color='red', s=10)  # s is the size of the dot
#     ax.axis('off')
#     plt.show()
#
#     return patches


def mortensen_single_frame(image_path, centroids_nm):

    # Load image - this should be an array of floats from 0-????
    # (let's just use mortensen's example data to set the max for now I guess)
    image = cv2.imread(image_path, 0)

    mortensen_example_data = np.loadtxt('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/mortensen/mortensen-example-data.txt')
    mort_min = np.min(mortensen_example_data)
    mort_max = np.max(mortensen_example_data)
    mort_range = mort_max - mort_min
    image = (image.astype(np.float64) / 255) * mort_range + mort_min

    # # Subtract the background roughly
    # image = np.clip(image - np.mean(image), 0, 255).astype(np.uint8)
    # image = image - np.mean(image)

    # Loop over all blobs
    x_list, y_list, phi_list, theta_list, covariance_list = [], [], [], [], []

    # for blob in blob_patches:
    for blob_index in range(len(centroids_nm)):

        # define patch of image to be considered
        # this is passed in the form of initpix and deltapix (centre pixel of patch, half width of patch)

        # # if using ground truth locations
        # x_nm = centroids_nm[blob_index][1]
        # y_nm = centroids_nm[blob_index][2]
        # if using slightly-wrong ground truth locations
        x_nm = centroids_nm[blob_index][1] + 100
        y_nm = centroids_nm[blob_index][2] + 100

        x_px = nm_to_px(x_nm, pixel_size_nm, image_size_px, 'x')
        y_px = nm_to_px(y_nm, pixel_size_nm, image_size_px, 'y')

        initpix = (int(x_px), int(y_px)) # Pixel coordinates of the center pixel in the area to be analyzed (patch should be centred on blob, so use that)

        image_patch = image[int(y_px - deltapix):int(y_px + deltapix), int(x_px - deltapix):int(x_px + deltapix)]

        image_with_patches = image.copy()

        # rescale to 0-255 for display with OpenCV
        image_with_patches = ((image_with_patches - np.min(image_with_patches)) /
                                    (np.max(image_with_patches) - np.min(image_with_patches)) * 255).astype(np.uint8)

        image_with_patches = cv2.cvtColor(image_with_patches, cv2.COLOR_GRAY2BGR)
        cv2.rectangle(image_with_patches, (int(x_px - deltapix), int(y_px - deltapix)), (int(x_px + deltapix), int(y_px + deltapix)), (255, 0, 0), 1)

        # # Plot patches
        # fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        # ax[0].imshow(image_with_patches, cmap='gray')
        # ax[0].axis('off')
        # ax[1].imshow(image_patch, cmap='gray')
        # ax[1].axis('off')
        # plt.show()

        start_blob = time.time()
        print(f"Analysing blob {blob_index+1}/{len(centroids_nm)}")

        # Create instance of MLEwT !!! this was inside the loop earlier !!!
        track = diPOLE_python3.MLEwT(peak_emission_wavelength,
                                     pixel_size_nm,
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

        x_relative_estimate_nm, y_relative_estimate_nm, theta_est, phi_est, cov_mat = track.Estimate(image)
        x_relative_estimates_nm.append(float(x_relative_estimate_nm))
        y_relative_estimates_nm.append(float(y_relative_estimate_nm))
        theta_list.append(float(theta_est))
        phi_list.append(float(phi_est))
        covariance_list.append(cov_mat)

        end_blob = time.time()
        elapsed_time_blob = end_blob - start_blob
        print(f"Time: {elapsed_time_blob:.4f} seconds on this blob")

    # # show on image
    # image_results = image_with_patches.copy()
    # image_size_nm = image_size_px * pixel_size_nm
    # # ground truth
    # true_x_nm = [centroid[1] for centroid in centroids_nm]
    # true_y_nm = [centroid[2] for centroid in centroids_nm]
    # # true_x_px = np.array([nm_to_px(x, pixel_size_nm, image.shape[1], 'x') for x in true_x_nm])
    # # true_y_px = np.array([nm_to_px(y, pixel_size_nm, image.shape[1], 'x') for y in true_y_nm]) # using 'x' because scatter uses bottom-left origin
    # # estimates
    # estimate_x_nm = true_x_nm[0] + x_relative_estimates_nm[0]
    # estimate_y_nm = true_y_nm[0] + y_relative_estimates_nm[0]
    # # estimate_x_px = np.array([nm_to_px(x, pixel_size_nm, image.shape[1], 'x') for x in estimate_x_nm])
    # # estimate_y_px = np.array([nm_to_px(y, pixel_size_nm, image.shape[1], 'x') for y in estimate_y_nm]) # using 'x' because scatter uses bottom-left origin
    # # print localisation error
    # print('ground truth')
    # for x, y in zip(true_x_nm, true_y_nm):
    #     print('(', x, y, ')')
    # print('estimates')
    # for x, y in zip(estimate_x_nm, estimate_y_nm):
    #     print('(', x, y, ')')
    # print('difference between ground truth and estimates')
    # for est_x, true_x, est_y, true_y in zip(estimate_x_nm, true_x_nm, estimate_y_nm, true_y_nm):
    #     print(np.abs(est_x - true_x), np.abs(est_y - true_y))
    # # plot
    # fig, ax = plt.subplots(figsize=(5, 5))
    # ax.imshow(image_results, extent=[-image_size_nm/2, image_size_nm/2, -image_size_nm/2, image_size_nm/2], cmap='gray')
    # ax.scatter(true_x_nm, true_y_nm, color='red', s=10)  # s is the size of the dot
    # ax.scatter(estimate_x_nm, estimate_y_nm, color='green', s=10)  # s is the size of the dot
    # ax.axis('off')
    # plt.show()

    return x_relative_estimates_nm, y_relative_estimates_nm, theta_list, phi_list, covariance_list


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
pixel_size_nm = 51.2 # Pixel width (nm per px)

# PSF parameters
photon_number = 5000.0 # Number of photons in image
background_level = 1.0 # Background level?
mu = 0.0001 # Initial guess location I think (nm)
nu = 0.0001 # Initial guess location I think (nm)
phi = 0.0001 # inclination
theta = 0.0001 # azimuthal angle
deltaz = 0.0 # Distance from design focal plane

# EMCCD parameters
inverse_gain = 1./100.
sigma_noise = 2 #12.6
Sfloor = 300.0
gain = 1.0 / inverse_gain

image_size_px = 101
patch_width_px = 11 # size of NxN pixel patch around blob centroid to consider

# --------------------

# Initial thunderstorm run to get blob location

# Load data
# frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/ASIL240923C05_50.ome.tif-frames/'
# results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/thunderstorm_results.csv'
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/mortensen/hinterer_sim/'
# results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/hinterer_results.csv'

# # find centroids using gaussian fitting thunderstorm
# centroids_image_coords = blob_detect_all_frames(frames_dir, results_path, pixel_size_nm)
# or for now give it the ground truth from some simulation. this is an array like [fram numbers array, x coords array, y coords array]

# # 9 spot image
# x_array_nm = [-1.8099*10**3, 0*10**3, 1.8099*10**3, -1.8099*10**3, 0*10**3, 1.8099*10**3, -1.8099*10**3, 0*10**3, 1.8099*10**3]
# y_array_nm = [1.8099*10**3, 1.8099*10**3, 1.8099*10**3, 0*10**3, 0*10**3, 0*10**3, -1.8099*10**3, -1.8099*10**3, -1.8099*10**3]

# 4 spot image
x_array_nm = [-1.2494*10**3]
y_array_nm = [-671.1673]

# # 1 spot image
# x_array_nm = [-1.8099*10**3]
# y_array_nm = [1.8099*10**3]

f_array = array = np.ones(len(x_array_nm))

centroids_nm = [(int(f), int(x), int(y)) for f, x, y in zip(f_array, x_array_nm, y_array_nm)]
# centroids_image_coords = [(int(f), int(nm_to_px(x, pixel_size_nm, image_size_px, 'x')), int(nm_to_px(y, pixel_size_nm, image_size_px, 'y'))) for f, x, y in zip(f_array, x_array_nm, y_array_nm)]

# get frames
frame_paths = [os.path.join(frames_dir, f) for f in os.listdir(frames_dir) if f.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'))]

# Initial guess params
initvals = np.array([mu, nu, background_level, photon_number, phi, theta, deltaz]) # initial PSF values
deltapix = patch_width_px / 2 # Half-width of the pixel area to be analyzed

# Mortensen run on each blob in each frame
x_relative_estimates_nm, y_relative_estimates_nm, theta_ests, phi_ests, covariance_ests = [], [], [], [], []
for i, frame_path in enumerate(frame_paths, 1):

    print("----------")
    print(f"FRAME {i}")
    start_frame = time.time()

    single_frame_results = list(mortensen_single_frame(frame_path, centroids_nm))

    x_relative_estimates_nm.append(single_frame_results[0])
    y_relative_estimates_nm.append(single_frame_results[1])
    theta_ests.append(single_frame_results[2])
    phi_ests.append(single_frame_results[3])
    # covariance_ests.append(single_frame_results[4])

    end_frame = time.time()
    elapsed_time_frame = end_frame - start_frame
    elapsed_time_frame = elapsed_time_frame/60
    print(f"Time: {elapsed_time_frame:.4f} minutes on this frame")

# make sure list is flat, because it needs to be for results table
# x_relative_estimates_nm = [item for sublist in x_relative_estimates_nm for item in sublist]
# y_relative_estimates_nm = [item for sublist in y_relative_estimates_nm for item in sublist]
# theta_ests = [item for sublist in theta_ests for item in sublist]
# phi_ests = [item for sublist in phi_ests for item in sublist]
# covariance_ests = [item for sublist in covariance_ests for item in sublist]

# # 9 spot ground truth
# ground_x = [-1.8099*10**3, 0, 1.8099*10**3, -1.8099*10**3, 0, 1.8099*10**3, -1.8099*10**3, 0, 1.8099*10**3]
# ground_y = [1.8099*10**3, 1.8099*10**3, 1.8099*10**3, 0, 0, 0, -1.8099*10**3, -1.8099*10**3, -1.8099*10**3]
# ground_theta = [0, 0.1745, 0.3491, 0.5236, 0.6981, 0.8727, 1.0472, 1.2217, 1.3963]
# ground_phi = [0, 0, 0, 0, 0, 0, 0, 0, 0]

# 1 spot ground truth
ground_x = [-1.2494*10**3]
ground_y = [-671.1673]
ground_theta = [0]
ground_phi = [0]

print('ground truth (x,y)')
print('(', ground_x[0], ',', ground_y[0], ')')
print('estimate (x,y)')
print('(', round(ground_x[0] + x_relative_estimates_nm[0],2), ',', round(ground_y[0] + y_relative_estimates_nm[0],2), ')')
print('error (x,y)')
print('(', round(x_relative_estimates_nm[0],2), ',', round(y_relative_estimates_nm[0],2), ')')

print('ground truth (θ,ϕ)')
print('(', ground_theta[0], ',', ground_phi[0], ')')
print('estimate (θ,ϕ)')

print('(', ground_theta[0] + round(theta_ests[0][0], 2), ',', ground_phi[0] + round(phi_ests[0][0],2), ')')
print('error (θ,ϕ)')
print('(', round(theta_ests[0][0],2), ',', round(phi_ests[0][0],2), ')')

# show on image
image_results = cv2.imread('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/mortensen/hinterer_sim/realistic-res.tif', 0)
image_results = cv2.cvtColor(image_results, cv2.COLOR_GRAY2RGB)

image_size_nm = image_size_px * pixel_size_nm

fig, ax = plt.subplots(figsize=(5, 5))
ax.imshow(image_results, extent=[-image_size_nm / 2, image_size_nm / 2, -image_size_nm / 2, image_size_nm / 2],
          cmap='gray')
ax.scatter(ground_x[0], ground_y[0], color='red', s=10)  # s is the size of the dot
ax.scatter(ground_x[0] + x_relative_estimates_nm[0], ground_y[0] + y_relative_estimates_nm[0], color='green', s=10)  # s is the size of the dot
ax.axis('off')
plt.show()

# print(len(covariance_ests))

# # generate a thunderstorm-style results table with new x,y localisations
# mortensen_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/mortensen_results.csv'
# # output_img_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/reconstruction.png'
#
# df = pd.read_csv(results_path)
# if len(x_relative_estimates_nm) != len(df) or len(y_relative_estimates_nm) != len(df):
#     raise ValueError("The length of the new x and y arrays must match the number of rows in the CSV file.")
# df['x [nm]'] += x_relative_estimates_nm
# df['y [nm]'] += y_relative_estimates_nm
# df.to_csv(mortensen_results_path, index=False)

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
