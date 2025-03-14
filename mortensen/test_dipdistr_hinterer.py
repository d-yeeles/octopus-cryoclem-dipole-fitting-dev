"""
Python implementation of the Hinterer PSF model.

This module provides a Python 3 compatible implementation of the PSF model originally
implemented in MATLAB in the 'hinterer' folder. It follows the same API as the 
Mortensen implementation to allow for interchangeable use.

Date: February 25, 2025
"""

import numpy as np
from numpy import *
import math
from scipy.special import jn

class Length:
    """Python implementation of the MATLAB Length class."""
    
    UNITS = {
        'nm': 1e-9,  # nanometer
        'mm': 1e-3,  # millimeter
        'm': 1.0      # meter (base unit)
    }
    
    def __init__(self, value, unit='nm'):
        """Initialize length value with unit."""
        self.value = np.array(value)
        self.unit = unit
        self._meter_value = self.value * self.UNITS[unit]
        
    def inMeter(self):
        """Return value in meters."""
        return self._meter_value
    
    def inNanometer(self):
        """Return value in nanometers."""
        return self._meter_value / self.UNITS['nm']
    
    def inPixels(self, pixel_size):
        """Return value in pixels given a pixel size."""
        return self._meter_value / pixel_size.inMeter()


class Dipole:
    """Python implementation of the MATLAB Dipole class."""
    
    def __init__(self, inclination=0, azimuth=0):
        """Initialize dipole with inclination and azimuth angles."""
        self.inclination = inclination
        self.azimuth = azimuth
        
    def getDipoleVector(self):
        """Return dipole orientation vector."""
        return [self.inclination, self.azimuth]


class Mask:
    """Python implementation of the MATLAB Mask class."""
    
    def __init__(self, n_grid=129, spacing=1, mode='FFT'):
        """Initialize mask with grid size and spacing."""
        self.n_grid = n_grid
        self.spacing = spacing
        self.mode = mode
        self.values = self._calculate_values()
        self.radius = self._calculate_radius()
        
    def _calculate_values(self):
        """Calculate mask values."""
        N = [self.n_grid, self.n_grid]
        s = [self.spacing, self.spacing]
        
        if self.mode == 'exact':
            x = np.linspace(-(N[0]-1)/2, (N[0]-1)/2, N[0]) * s[0]
            y = np.linspace(-(N[1]-1)/2, (N[1]-1)/2, N[1]) * s[1]
        elif self.mode == 'FFT':
            x = np.arange(np.floor(-N[0]/2 + 0.5), np.floor(N[0]/2 - 0.5) + 1) * s[0]
            y = np.arange(np.floor(-N[1]/2 + 0.5), np.floor(N[1]/2 - 0.5) + 1) * s[1]
            
        X, Y = np.meshgrid(y, x)
        radius = np.sqrt(X**2 + Y**2)
        return (radius**2 <= (min(N[0]/2*s[0], N[1]/2*s[1]))**2).astype(float)
    
    def _calculate_radius(self):
        """Calculate radius values."""
        x = np.linspace(-1, 1, self.n_grid)
        X, Y = np.meshgrid(x, x)
        return np.sqrt(X**2 + Y**2)
    
    def getPolarCoordinates(self):
        """Return polar coordinates of the mask."""
        x = np.linspace(-1, 1, self.n_grid)
        X, Y = np.meshgrid(x, x)
        angle, radius = np.flipud(np.arctan2(Y, X) + np.pi), np.sqrt(X**2 + Y**2)
        return radius, angle


class BackFocalPlane:
    """Python implementation of the MATLAB BackFocalPlane class."""
    
    def __init__(self, psf_obj):
        """Initialize back focal plane with PSF object."""
        self.n_grid = psf_obj.n_discretization_bfp
        self.unit_k_space = psf_obj.unit_k_space
        self.electric_field = self._calculate_electric_field(psf_obj)
        
    def _calculate_electric_field(self, psf_obj):
        """Calculate electric field in the back focal plane."""
        # Extract parameters
        RI = psf_obj.refractive_indices
        dipole = psf_obj.dipole
        h_intermediate = psf_obj.height_intermediate_layer.inMeter()
        pos = psf_obj.position.inMeter()
        z = pos[2]
        focal_length = psf_obj.objective_focal_length.inMeter()
        mu = 1e-12  # Dipole moment magnitude
        
        # Pre-calculations
        if len(RI) == 1:
            RI = [RI[0], RI[0], RI[0]]
            
        # Coordinates in the objective pupil
        pupil_mask = Mask(self.n_grid, self.unit_k_space)
        _, mask_angle = pupil_mask.getPolarCoordinates()
        PHI3 = np.fliplr(mask_angle) - np.pi
        
        # Wavenumbers
        k0 = 2 * np.pi / psf_obj.wavelength.inMeter()  # vacuum
        k1 = k0 * RI[0]  # medium 1 (typically water)
        k2 = k0 * RI[1]  # medium 2 (intermediate layer)
        k3 = k0 * RI[2]  # medium 3 (immersion medium)
        
        # Angles in different media
        Kr = pupil_mask.radius
        THETA1 = np.arccos(np.sqrt(1 - (RI[2]/RI[0] * Kr/k3)**2)) * pupil_mask.values
        THETA2 = np.arccos(np.sqrt(1 - (RI[2]/RI[1] * Kr/k3)**2)) * pupil_mask.values
        THETA3 = np.arcsin(Kr/k3) * pupil_mask.values
        
        # Calculations according to Axelrod, 2012
        # Cosines of angles
        CTHETA1 = np.cos(THETA1)
        CTHETA2 = np.cos(THETA2)
        CTHETA3 = np.cos(THETA3)
        
        # Fresnel coefficients
        tp12 = 2*RI[0]*CTHETA1/(RI[0]*CTHETA2+RI[1]*CTHETA1)
        tp23 = 2*RI[1]*CTHETA2/(RI[1]*CTHETA3+RI[2]*CTHETA2)
        
        ts12 = 2*RI[0]*CTHETA1/(RI[0]*CTHETA1+RI[1]*CTHETA2)
        ts23 = 2*RI[1]*CTHETA2/(RI[1]*CTHETA2+RI[2]*CTHETA3)
        
        rp12 = (RI[1]*CTHETA1-RI[0]*CTHETA2)/(RI[0]*CTHETA2+RI[1]*CTHETA1)
        rp23 = (RI[2]*CTHETA2-RI[1]*CTHETA3)/(RI[1]*CTHETA3+RI[2]*CTHETA2)
        
        rs12 = (RI[0]*CTHETA1-RI[1]*CTHETA2)/(RI[0]*CTHETA1+RI[1]*CTHETA2)
        rs23 = (RI[1]*CTHETA2-RI[2]*CTHETA3)/(RI[1]*CTHETA2+RI[2]*CTHETA3)
        
        # Fresnel coefficients for three-layer system
        tp = tp12 * tp23 * np.exp(1j*k2*h_intermediate*CTHETA2) / (1 + rp12 * rp23 * np.exp(2j*k2*h_intermediate*CTHETA2))
        ts = ts12 * ts23 * np.exp(1j*k2*h_intermediate*CTHETA2) / (1 + rs12 * rs23 * np.exp(2j*k2*h_intermediate*CTHETA2))
        
        # Dipole projections
        mu_p = mu * np.sin(dipole.inclination) * np.cos(dipole.azimuth - PHI3)
        mu_s = mu * np.sin(dipole.inclination) * np.sin(dipole.azimuth - PHI3)
        mu_z = mu * np.cos(dipole.inclination)
        
        # Prefactor C
        C = (k3**2 * np.exp(1j*k3*focal_length) * CTHETA3) / (focal_length * RI[0]) * \
            np.exp(-1j*k3*h_intermediate*CTHETA3) * \
            np.exp(1j*k1*CTHETA1*z)
            
        # Electric field components in layer 3
        E3p = C * tp * CTHETA3 * (mu_p/RI[2] + mu_z*np.sin(THETA3)/CTHETA1)
        E3s = C * ts * (mu_s/(RI[2]/CTHETA1))
        E3z = C * tp * np.sin(THETA3) * (mu_p/RI[2] + mu_z*np.sin(THETA3)/CTHETA1)
        
        # Influence of objective
        apodization_factor = 1 / np.sqrt(CTHETA3) * pupil_mask.values
        E_BFP_p = (E3p*CTHETA3 + E3z*np.sin(THETA3)) * apodization_factor
        E_BFP_s = E3s * apodization_factor
        
        # Coordinate transformation into x-and y-polarization
        E_BFP_x = np.cos(PHI3)*E_BFP_p - np.sin(PHI3)*E_BFP_s
        E_BFP_y = np.sin(PHI3)*E_BFP_p + np.cos(PHI3)*E_BFP_s
        
        return {'x': E_BFP_x, 'y': E_BFP_y}
        
class BackFocalPlane_gaussian:
    """Python implementation of the MATLAB BackFocalPlane_gaussian class."""
    
    def __init__(self, psf_obj):
        """Initialize gaussian back focal plane with PSF object."""
        self.n_grid = psf_obj.n_discretization_bfp
        self.unit_k_space = psf_obj.unit_k_space
        self.electric_field = self._calculate_electric_field(psf_obj)
        
    def _calculate_electric_field(self, psf_obj):
        """Calculate gaussian electric field in the back focal plane."""
        # Coordinates in the objective pupil
        pupil_mask = Mask(self.n_grid, self.unit_k_space)
        Kr, _ = pupil_mask.getPolarCoordinates()
        
        # Define Gaussian parameters
        sigma = 0.3  # Standard deviation of the Gaussian
        amplitude = 1  # Peak amplitude of the Gaussian
        
        # Gaussian Electric Field
        gaussian_field = amplitude * np.exp(-Kr**2 / (2 * sigma**2))
        
        # Output Gaussian electric field for both x and y components
        return {'x': gaussian_field, 'y': gaussian_field}


class ChirpZTransform:
    """Python implementation of the MATLAB ChirpZTransform class."""
    
    def __init__(self, psf_obj):
        """Initialize chirp-z transform with PSF object."""
        # Get dimensions of image and k-space
        n_pixels_image = psf_obj.n_pixels * psf_obj.oversampling
        n_pixels_bfp = psf_obj.n_discretization_bfp
        self.n_pixels_padded = n_pixels_image + n_pixels_bfp - 1
        
        # Create grid
        x = np.arange(np.floor(-self.n_pixels_padded/2 + 0.5), 
                     np.floor(self.n_pixels_padded/2 - 0.5) + 1)
        y = np.arange(np.floor(-self.n_pixels_padded/2 + 0.5), 
                     np.floor(self.n_pixels_padded/2 - 0.5) + 1)
        
        # Create kernel (quadratic convolution phase kernel)
        alpha = psf_obj.unit_object_space * psf_obj.unit_k_space / (2*np.pi)
        X, Y = np.meshgrid(y, x)
        self.kernel = np.exp(-1j*alpha*np.pi*(X**2 + Y**2))
        
        # Fourier transform of kernel
        self.fourier_kernel = np.fft.fft2(np.fft.ifftshift(self.kernel))
        
    def apply(self, psf_obj, E_in):
        """Apply chirp-z transform to input electric field."""
        # Complex conjugate of kernel
        conj_kernel = np.conjugate(self.kernel)
        
        # Electric field times complex conjugate of kernel
        f = self._embed_array_2d(E_in, self.n_pixels_padded) * conj_kernel
        
        # Fourier transform of f
        F = np.fft.fft2(np.fft.ifftshift(f))
        
        # Convolving f with kernel (inverse fft of multiplication)
        convolution = np.fft.fftshift(np.fft.ifft2(F * self.fourier_kernel))
        
        # Final multiplication by factor
        E_out = self.n_pixels_padded * conj_kernel * convolution
        
        # Crop to image size
        E_out = self._crop_array_2d(E_out, psf_obj.n_pixels * psf_obj.oversampling)
        
        return E_out
    
    def _embed_array_2d(self, array_in, target_size, pad_val=0):
        """Embed a 2D array inside a bigger canvas."""
        if isinstance(target_size, (int, float)):
            target_size = [target_size, target_size]
            
        size_in = array_in.shape
        size_diff = np.array(target_size) - np.array(size_in)
        
        if np.any(size_diff < 0):
            raise ValueError("Target size for embedding must be greater than size of input array!")
        
        pad_above = int(np.ceil(size_diff[0]/2))
        pad_below = int(np.floor(size_diff[0]/2))
        pad_left = int(np.ceil(size_diff[1]/2))
        pad_right = int(np.floor(size_diff[1]/2))
        
        return np.pad(array_in, ((pad_above, pad_below), (pad_left, pad_right)), 
                     mode='constant', constant_values=pad_val)
    
    def _crop_array_2d(self, array_in, target_size):
        """Crop a 2D array to target size."""
        if isinstance(target_size, (int, float)):
            target_size = [target_size, target_size]
            
        size_in = array_in.shape
        size_diff = np.array(size_in) - np.array(target_size)
        
        if np.any(size_diff < 0):
            raise ValueError("Target size for cropping must be smaller than size of input array!")
        
        start_row = int(np.floor(size_diff[0]/2))
        end_row = int(start_row + target_size[0])
        start_col = int(np.floor(size_diff[1]/2))
        end_col = int(start_col + target_size[1])
        
        return array_in[start_row:end_row, start_col:end_col]


class dipdistr_hinterer:
    """
    Python implementation of the Hinterer PSF model.
    
    This class calculates the theoretical point spread function (PSF) for fixed dipoles
    following the implementation in the Hinterer MATLAB codebase. It follows the same API
    as the Mortensen implementation to allow for interchangeable use.
    
    Parameters:
    -----------
    wl : float
        Wavelength [nm] of emission peak in buffer.
    NA : float
        Numerical aperture of the objective.
    n : float
        Refractive index of glass/immersion oil.
    n0 : float
        Refractive index of buffer (water).
    M : float
        Magnification of the objective.
        
    Methods:
    --------
    PSF_exact : float
        Takes the coordinates in the image, the angles of the probe, and the distance 
        from the focus as arguments. Returns the value of the exact PSF at that point.
        
    PSF_approx : float
        Takes the coordinates in the image, the angles of the probe, and the distance 
        from the focus as arguments. Returns the value of the approximated PSF at that point.
    """
    
    def __init__(self, wl, n, n0, M, NA):
        """Initialize PSF model with given parameters."""
        self.wl = wl            # Wavelength in nm (emission peak in buffer)
        self.NA = NA            # Numerical aperture
        self.M = M              # Magnification
        self.n = n              # Refractive index of glass/oil
        self.n0 = n0            # Refractive index of sample medium (water)
        
        # Calculate normalized PSF
        self._setup_psf_parameters()
        
        # Store position and focus
        self.rho = 0.0
        self.deltaz = 0.0
    
    def _setup_psf_parameters(self):
        """Set up PSF parameters."""
        # Create PSF object
        self.psf_params = {
            'n_pixels': 17,
            'dipole': Dipole(0, 0),
            'position': Length([0, 0, 0], 'nm'),
            'n_photons': 1e5,
            'wavelength': Length(self.wl, 'nm'),
            'defocus': Length(0, 'nm'),
            'astigmatism': 0,
            'objective_NA': self.NA,
            'objective_focal_length': Length(3, 'mm'),
            'refractive_indices': [self.n0, self.n, self.n],
            'height_intermediate_layer': Length(0, 'mm'),
            'n_discretization_bfp': 129,
            'pixel_size': Length(108, 'nm'),
            'pixel_sensitivity_mask': np.ones((9, 9)),
            'background_noise': 0,
            'oversampling': 9,
            'unit_object_space': 108e-9 / 9,  # pixelSize/oversampling
            'unit_k_space': 1
        }
        
        # Initialize back focal plane calculations
        self.bfp = None
        self.bfp_gaussian = None
        self.chirp_z_transform = None
        self.defocus = None
        self.pupil_mask = None
        
def SetFocus(self, deltaz):
        """Set the focus distance."""
        self.deltaz = deltaz
        
    def GetFocus(self):
        """Get the focus distance."""
        return self.deltaz
    
    def SetRho(self, rho):
        """Set the radial position."""
        self.rho = rho
        
    def GetRho(self):
        """Get the radial position."""
        return self.rho
    
    def _calculate_psf(self, x, y, alpha, beta, deltaz, use_gaussian=False):
        """
        Calculate the PSF for given coordinates and angles.
        
        Parameters:
        -----------
        x, y : float
            Coordinates in the image plane.
        alpha, beta : float
            Angles of the dipole orientation.
        deltaz : float
            Distance from focus.
        use_gaussian : bool
            Whether to use Gaussian approximation.
            
        Returns:
        --------
        psf_value : float
            The PSF value at the given point.
        """
        # Update PSF parameters
        self.psf_params['dipole'] = Dipole(beta, alpha + np.pi)  # Adjust orientation
        self.psf_params['position'] = Length([x, y, 0], 'nm')
        self.psf_params['defocus'] = Length(deltaz, 'nm')
        
        # Initialize components if needed
        if self.pupil_mask is None:
            self.pupil_mask = Mask(self.psf_params['n_discretization_bfp'])
        
        # Calculate electric field in back focal plane
        if use_gaussian:
            if self.bfp_gaussian is None:
                self.bfp_gaussian = BackFocalPlane_gaussian(self.psf_params)
            E_bfp = self.bfp_gaussian.electric_field
        else:
            if self.bfp is None:
                self.bfp = BackFocalPlane(self.psf_params)
            E_bfp = self.bfp.electric_field
        
        # Apply chirp-z transform to get image plane field
        if self.chirp_z_transform is None:
            self.chirp_z_transform = ChirpZTransform(self.psf_params)
            
        # Calculate intensities
        Ix = np.abs(self.chirp_z_transform.apply(self.psf_params, E_bfp['x']))**2
        Iy = np.abs(self.chirp_z_transform.apply(self.psf_params, E_bfp['y']))**2
        
        # Sum intensities
        psf = Ix + Iy
        
        # Normalize and apply photon count
        total_intensity = np.sum(psf)
        psf = psf / total_intensity * self.psf_params['n_photons']
        
        # Return intensity at center pixel
        center = psf.shape[0] // 2
        return psf[center, center]
    
    def PSF_exact(self, x, y, alpha, beta, deltaz):
        """
        Calculate the exact PSF value for given parameters.
        
        Parameters:
        -----------
        x, y : float
            Coordinates in the image plane.
        alpha, beta : float
            Angles of the dipole orientation.
        deltaz : float
            Distance from focus.
            
        Returns:
        --------
        psf_value : float
            The PSF value at the given point.
        """
        return self._calculate_psf(x, y, alpha, beta, deltaz, use_gaussian=False)
    
    def PSF_approx(self, x, y, alpha, beta, deltaz):
        """
        Calculate the approximate PSF value for given parameters.
        
        Parameters:
        -----------
        x, y : float
            Coordinates in the image plane.
        alpha, beta : float
            Angles of the dipole orientation.
        deltaz : float
            Distance from focus.
            
        Returns:
        --------
        psf_value : float
            The PSF value at the given point.
        """
        return self._calculate_psf(x, y, alpha, beta, deltaz, use_gaussian=True)
        

