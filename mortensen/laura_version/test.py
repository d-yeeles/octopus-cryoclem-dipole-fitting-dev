import diPOLE_mixed
from pylab import *
import numpy as np
import cv2
import os
import time
import tifffile

def nm_to_px(coord_nm, pixel_size_nm, image_size_px, x_or_y_flag):
    if x_or_y_flag == 'x':
        coord_px = (coord_nm / pixel_size_nm + (image_size_px + 1) / 2) - 1
    elif x_or_y_flag == 'y':
        coord_px = image_size_px - (coord_nm / pixel_size_nm + image_size_px / 2) - 1
    return coord_px

def mortensen_single_frame(image_path):

    # Load image - this should be an array of floats from 0-????
    # (let's just use mortensen's example data to set the max for now I guess)
    # image = cv2.imread(image_path, 0)
    image = tifffile.imread(image_path)
    mortensen_example_data = np.loadtxt('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/mortensen/mortensen-example-data.txt')
    mort_min = np.min(mortensen_example_data)
    mort_max = np.max(mortensen_example_data)
    mort_range = mort_max - mort_min
    # image = (image.astype(np.float64) / 255) * mort_range + mort_min

    # fig, ax = plt.subplots(figsize=(5, 5))
    # ax.imshow(image, extent=[-image_size_nm / 2, image_size_nm / 2, -image_size_nm / 2, image_size_nm / 2],
    #           cmap='gray')
    # ax.axis('off')
    # plt.show()

    # Loop over all blobs
    x_list, y_list, phi_list, theta_list, photon_list = [], [], [], [], []

    start_blob = time.time()

    # Create instance of MLEwT
    track = diPOLE_mixed.MLEwT(peak_emission_wavelength,
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

    phi_est, theta_est, x_est, y_est, photon_est = track.Estimate(image)
    x_list.append(float(x_est))
    y_list.append(float(y_est))
    theta_list.append(float(theta_est))
    phi_list.append(float(phi_est))
    photon_list.append(float(photon_est))

    end_blob = time.time()
    elapsed_time_blob = end_blob - start_blob
    print(f"Time: {elapsed_time_blob:.2f} seconds on this blob")

    return x_list[0], y_list[0], theta_list[0], phi_list[0], photon_list[0]



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
photon_number = 1e6 # Number of photons in image
background_level = 1.0 # Background level
mu = 12.345 # Initial guess location I think (nm)
nu = 12.345 # Initial guess location I think (nm)
phi = 0. # inclination
theta = 0. # azimuthal angle
deltaz = 0. # Distance from design focal plane

# EMCCD parameters
inverse_gain = 1./100.
sigma_noise = 2 #12.6
Sfloor = 300.0
gain = 1.0 / inverse_gain

image_size_px = 19
image_size_nm = image_size_px * pixel_size_nm
patch_width_px = image_size_px # size of NxN pixel patch around blob centroid to consider

# --------------------

# Load data
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/mortensen/hinterer_sim/'

# get frames
frame_paths = [os.path.join(frames_dir, f) for f in os.listdir(frames_dir) if f.endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'))]

# Initial guess params
initvals = np.array([mu, nu, background_level, photon_number, phi, theta, deltaz]) # initial PSF values
deltapix = round(image_size_px/2) #patch_width_px / 2 # Half-width of the pixel area to be analyzed
initpix = (deltapix, deltapix)

# Mortensen run on each blob in each frame
x_ests, y_ests, theta_ests, phi_ests, photon_ests = [], [], [], [], []

for i, frame_path in enumerate(frame_paths, 1):

    print("----------")
    print(f"FRAME {i}")
    start_frame = time.time()

    single_frame_results = list(mortensen_single_frame(frame_path))

    x_ests.append(single_frame_results[0])
    y_ests.append(single_frame_results[1])
    theta_ests.append(single_frame_results[2])
    phi_ests.append(single_frame_results[3])
    photon_ests.append(single_frame_results[4])

    end_frame = time.time()
    elapsed_time_frame = end_frame - start_frame
    elapsed_time_frame = elapsed_time_frame/60
    print(f"Time: {elapsed_time_frame:.2f} minutes on this frame")

x_true = [-1.485782e+01]
y_true = [-2.518784e+01]
theta_true = [0]
phi_true = [0]

print('estimates:')
print(x_ests, y_ests, theta_ests, phi_ests)

print('errors:')
print([a - b for a, b in zip(x_true, x_ests)])
print([a - b for a, b in zip(y_true, y_ests)])
print([(a - b)*180/np.pi for a, b in zip(theta_true, theta_ests)])
print([(a - b)*180/np.pi for a, b in zip(phi_true, phi_ests)])

# show on image
image_results = tifffile.imread('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/mortensen/hinterer_sim/sim_inc000_az000_run1.tif')
# image_results = cv2.cvtColor(image_results, cv2.COLOR_GRAY2RGB)

fig, ax = plt.subplots(figsize=(5, 5))
ax.imshow(image_results, extent=[-image_size_nm / 2, image_size_nm / 2, -image_size_nm / 2, image_size_nm / 2],
          cmap='gray')
ax.scatter(x_true[0], y_true[0], color='red', s=10)  # s is the size of the dot
ax.scatter(x_ests[0], y_ests[0], color='yellow', s=10)
ax.scatter(x_ests_reparam[0], y_ests_reparam[0], color='green', s=10)
ax.axis('off')
plt.show()
