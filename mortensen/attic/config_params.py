import numpy

# Experimental parameters
wl = 580.0 # Peak emission wavelength
n = 1.52 # RI of immersion oil/optics
n0 = 1.33 # RI of buffer
NA = 1.49 # Numerical aperture of objective
M = 250.0 # Magnification of composite microscope
deltax = 23.4 # Pixel width (nm per px)

# PSF parameters
N = 10000.0 # Photon number?
b = 5.0 # Background level?
mu = 0.1 # Probe location?
nu = 0.1 # ...
phi = 2 * numpy.pi / 3.0 # inclination
theta = 0.5 # azimuthal angle
deltaz = -30.0 # Distance from design focal plane

# EMCCD parameters
inversegain = 0.0089
sigmanoise = 12.6
Sfloor = 300.0
gain = 1.0 / inversegain

npix = 12 # size of NxN pixel patch around blob centroid to consider
