import numpy as np
from scipy.special import jn, gamma, erfc, i1
from scipy.optimize import fmin_powell
import time
import matplotlib.pyplot as plt


class Dipole:
    """Represents a dipole with orientation angles"""

    def __init__(self, inclination=0, azimuth=0):
        self.inclination = inclination  # Polar angle
        self.azimuth = azimuth  # Azimuthal angle


class Mask:
    """Recreates the Mask class from MATLAB"""

    def __init__(self, n_grid, spacing=1):
        self.n_grid = n_grid
        self.spacing = spacing
        self.values = np.zeros((n_grid, n_grid))

        # Calculate mask values
        x = np.linspace(-1, 1, n_grid)
        X, Y = np.meshgrid(x, x)
        self.X = X  # Store for later use
        self.Y = Y  # Store for later use
        self.radius = np.sqrt(X ** 2 + Y ** 2)
        self.values = (self.radius <= 1.0).astype(float)

    def get_polar_coordinates(self):
        """Return polar coordinates of mask points"""
        # Fix: Correctly calculate angles from X and Y coordinates
        angle = np.arctan2(self.Y, self.X)
        return angle, self.radius


class ZernikePolynomials:
    """Implements Zernike polynomials for aberration modeling"""

    @staticmethod
    def get_instance(mask):
        """Singleton-like instance getter"""
        return ZernikePolynomials(mask)

    def __init__(self, mask):
        self.mask = mask
        self.polynomials = {}

    def get_aberration(self, indices, coefficients):
        """Calculate weighted sum of Zernike polynomials"""
        aberration = np.zeros_like(self.mask.values)

        for idx, coeff in zip(indices, coefficients):
            if idx not in self.polynomials:
                self.polynomials[idx] = self._calculate_polynomial(idx)
            aberration += coeff * self.polynomials[idx]

        return aberration

    def _calculate_polynomial(self, idx):
        """Calculate a specific Zernike polynomial using the Noll index"""
        n, m = self._noll_to_zernike_indices(idx)
        return self.zernike(n, m, self.mask.radius, np.arctan2(self.mask.Y, self.mask.X)) * self.mask.values

    def _noll_to_zernike_indices(self, j):
        """
        Convert Noll index to (n,m) Zernike indices
        j: Noll index (integer starting at 1)
        Returns: tuple of (n,m) indices where:
          n: radial (non-negative integer)
          m: azimuthal (positive or negative integer with |m| <= n and n-|m| even)
        """
        if j <= 0:
            raise ValueError("Noll indices start at 1")

        n = 0
        j1 = j - 1
        while j1 > n:
            n += 1
            j1 -= n

        m = (-1) ** (j % 2) * ((n % 2) + 2 * int((j1 + ((n + 1) % 2)) / 2))

        return n, m

    def zernike(self, n, m, rho, theta):
        """
        Calculate Zernike polynomial Z_n^m(rho, theta)
        n: radial index (non-negative integer)
        m: azimuthal index (integer with |m| <= n and n-|m| even)
        rho: radial coordinate (normalized to [0,1])
        theta: azimuthal coordinate in radians
        Returns: Zernike polynomial values at coordinates (rho, theta)
        """
        m_abs = abs(m)

        # Validate inputs
        if n < 0:
            raise ValueError("Radial index n must be non-negative")
        if m_abs > n:
            raise ValueError("Azimuthal index |m| must be <= n")
        if (n - m_abs) % 2 != 0:
            raise ValueError("n - |m| must be even")

        # Calculate radial component
        R = self._radial_polynomial(n, m_abs, rho)

        # Calculate angular component
        if m > 0:
            A = np.cos(m * theta)
        elif m < 0:
            A = np.sin(m_abs * theta)
        else:  # m == 0
            A = np.ones_like(theta)

        # Normalization factor
        if m == 0:
            norm = np.sqrt(n + 1)
        else:
            norm = np.sqrt(2 * (n + 1))

        return R * A / norm

    def _radial_polynomial(self, n, m, rho):
        """
        Calculate radial polynomial R_n^m(rho)
        n: radial index (non-negative integer)
        m: absolute value of azimuthal index (|m| <= n and n-|m| even)
        rho: radial coordinate (normalized to [0,1])
        Returns: radial polynomial values at coordinates rho
        """
        if n == m:
            return rho ** n

        R = np.zeros_like(rho)

        for k in range((n - m) // 2 + 1):
            # Fix: Use scipy.special.factorial instead of recursive implementation
            coef = (-1) ** k * factorial(n - k) / (
                    factorial(k) * factorial((n + m) // 2 - k) * factorial((n - m) // 2 - k)
            )
            R += coef * rho ** (n - 2 * k)

        return R

    # Optional: Named Zernike polynomials for common aberrations
    @staticmethod
    def get_polynomial_name(idx):
        """Return common name for Zernike polynomial by Noll index"""
        names = {
            1: "Piston",
            2: "Tip (X-Tilt)",
            3: "Tilt (Y-Tilt)",
            4: "Defocus",
            5: "Oblique Astigmatism (45°)",
            6: "Vertical Astigmatism (0°)",
            7: "Vertical Coma",
            8: "Horizontal Coma",
            9: "Vertical Trefoil",
            10: "Oblique Trefoil",
            11: "Primary Spherical",
            12: "Vertical Secondary Astigmatism",
            13: "Oblique Secondary Astigmatism",
            14: "Vertical Quadrafoil",
            15: "Oblique Quadrafoil",
            16: "Secondary Vertical Coma",
            17: "Secondary Horizontal Coma",
            18: "Secondary Vertical Trefoil",
            19: "Secondary Oblique Trefoil",
            20: "Secondary Spherical",
            21: "Tertiary Vertical Astigmatism",
            22: "Tertiary Oblique Astigmatism"
        }

        return names.get(idx, f"Zernike polynomial {idx}")


class BackFocalPlane:
    """Calculates the back focal plane fields"""

    def __init__(self, psf_obj):
        self.electric_field = self.calculate_electric_field(psf_obj)

    def calculate_electric_field(self, psf_obj):
        """Calculate electric field in the back focal plane"""
        # Get parameters from PSF object
        RI = psf_obj.refractive_indices
        dipole = psf_obj.dipole
        h_intermediate = psf_obj.height_intermediate_layer
        pos = psf_obj.position
        z = pos[2]
        focal_length = psf_obj.objective_focal_length
        mu = 1e-12  # Dipole moment scale

        # Ensure RI has 3 values
        if len(RI) == 1:
            RI = [RI, RI, RI]

        # Create mask for pupil coordinates
        n_grid = psf_obj.n_discretization_bfp
        pupil_mask = Mask(n_grid)

        # Get polar coordinates from mask
        mask_angle, mask_radius = pupil_mask.get_polar_coordinates()

        # Fix: Don't use fliplr on the angles
        PHI3 = mask_angle

        # Wavenumbers in different media
        k0 = 2 * np.pi / psf_obj.wavelength  # vacuum
        k1 = k0 * RI[0]  # specimen
        k2 = k0 * RI[1]  # intermediate layer
        k3 = k0 * RI[2]  # immersion medium

        # Angles in different media
        Kr = pupil_mask.radius

        # Fix: Add safety checks to avoid numerical issues
        eps = 1e-10  # small value to avoid numerical issues
        THETA1 = np.arccos(np.sqrt(np.maximum(eps, 1 - (RI[2] / RI[0] * Kr / k3) ** 2))) * pupil_mask.values
        THETA2 = np.arccos(np.sqrt(np.maximum(eps, 1 - (RI[2] / RI[1] * Kr / k3) ** 2))) * pupil_mask.values
        THETA3 = np.arcsin(np.minimum(1 - eps, Kr / k3)) * pupil_mask.values

        # Cosines of angles
        CTHETA1 = np.cos(THETA1)
        CTHETA2 = np.cos(THETA2)
        CTHETA3 = np.cos(THETA3)

        # Fresnel coefficients
        tp12 = 2 * RI[0] * CTHETA1 / (RI[0] * CTHETA2 + RI[1] * CTHETA1)
        tp23 = 2 * RI[1] * CTHETA2 / (RI[1] * CTHETA3 + RI[2] * CTHETA2)

        ts12 = 2 * RI[0] * CTHETA1 / (RI[0] * CTHETA1 + RI[1] * CTHETA2)
        ts23 = 2 * RI[1] * CTHETA2 / (RI[1] * CTHETA2 + RI[2] * CTHETA3)

        rp12 = (RI[1] * CTHETA1 - RI[0] * CTHETA2) / (RI[0] * CTHETA2 + RI[1] * CTHETA1)
        rp23 = (RI[2] * CTHETA2 - RI[1] * CTHETA3) / (RI[1] * CTHETA3 + RI[2] * CTHETA2)

        rs12 = (RI[0] * CTHETA1 - RI[1] * CTHETA2) / (RI[0] * CTHETA1 + RI[1] * CTHETA2)
        rs23 = (RI[1] * CTHETA2 - RI[2] * CTHETA3) / (RI[1] * CTHETA2 + RI[2] * CTHETA3)

        # Fresnel coefficients for three-layer system
        tp = tp12 * tp23 * np.exp(1j * k2 * h_intermediate * CTHETA2) / (
                1 + rp12 * rp23 * np.exp(2j * k2 * h_intermediate * CTHETA2))
        ts = ts12 * ts23 * np.exp(1j * k2 * h_intermediate * CTHETA2) / (
                1 + rs12 * rs23 * np.exp(2j * k2 * h_intermediate * CTHETA2))

        # Dipole projections onto directions p, s and z
        mu_p = mu * np.sin(dipole.inclination) * np.cos(dipole.azimuth - PHI3)
        mu_s = mu * np.sin(dipole.inclination) * np.sin(dipole.azimuth - PHI3)
        mu_z = mu * np.cos(dipole.inclination)

        # Prefactor C
        C = (k3 ** 2 * np.exp(1j * k3 * focal_length) * CTHETA3) / (focal_length * RI[0]) \
            * np.exp(-1j * k3 * h_intermediate * CTHETA3) \
            * np.exp(1j * k1 * CTHETA1 * z)

        # Electric field components in layer 3
        E3p = C * tp * CTHETA3 * (mu_p / RI[2] + mu_z * np.sin(THETA3) / CTHETA1)
        E3s = C * ts * (mu_s / (RI[2] / CTHETA1))
        E3z = C * tp * np.sin(THETA3) * (mu_p / RI[2] + mu_z * np.sin(THETA3) / CTHETA1)

        # Influence of objective
        apodization_factor = 1 / np.sqrt(CTHETA3) * pupil_mask.values
        E_BFP_p = (E3p * CTHETA3 + E3z * np.sin(THETA3)) * apodization_factor
        E_BFP_s = E3s * apodization_factor

        # Coordinate transformation into x-and y-polarization
        E_BFP = {}
        E_BFP['x'] = np.cos(PHI3) * E_BFP_p - np.sin(PHI3) * E_BFP_s
        E_BFP['y'] = np.sin(PHI3) * E_BFP_p + np.cos(PHI3) * E_BFP_s

        return E_BFP


class ChirpZTransform:
    """Implements the Chirp-Z transform for field propagation"""

    def __init__(self, psf_obj):
        # Get dimensions
        n_pixels_image = psf_obj.n_pixels * psf_obj.oversampling
        n_pixels_bfp = psf_obj.n_discretization_bfp
        self.n_pixels_padded = n_pixels_image + n_pixels_bfp - 1

        # Create grid - make sure it has the correct dimensions
        n_pad = self.n_pixels_padded
        x = np.arange(-n_pad // 2, n_pad // 2) if n_pad % 2 == 0 else np.arange(-n_pad // 2, n_pad // 2 + 1)
        y = x.copy()  # Use the same range for both dimensions
        X, Y = np.meshgrid(x, y)

        # Create kernel
        alpha = psf_obj.unit_object_space * psf_obj.unit_k_space / (2 * np.pi)
        self.kernel = np.exp(-1j * alpha * np.pi * (X ** 2 + Y ** 2))

        # Fourier transform of kernel
        self.fourier_kernel = np.fft.fft2(np.fft.ifftshift(self.kernel))

    def apply(self, psf_obj, E_in):
        # Complex conjugate of kernel
        conj_kernel = np.conj(self.kernel)

        # Make sure the embedded array has the same dimensions as the kernel
        f = self._embed_array_2D(E_in, self.kernel.shape[0]) * conj_kernel

        # Fourier transform
        F = np.fft.fft2(np.fft.ifftshift(f))

        # Convolution (inverse fft of multiplication)
        convolution = np.fft.fftshift(np.fft.ifft2(F * self.fourier_kernel))

        # Final multiplication
        E_out = self.kernel.shape[0] * conj_kernel * convolution

        # Crop to image size
        E_out = self._crop_array_2D(E_out, psf_obj.n_pixels * psf_obj.oversampling)

        return E_out

    def _embed_array_2D(self, array, new_size):
        """Embed array in center of larger array"""
        old_size = array.shape

        # Create new array
        new_array = np.zeros((new_size, new_size), dtype=array.dtype)

        # Calculate offsets
        offset_x = (new_size - old_size[0]) // 2
        offset_y = (new_size - old_size[1]) // 2

        # Insert array into center
        new_array[offset_x:offset_x + old_size[0], offset_y:offset_y + old_size[1]] = array

        return new_array

    def _crop_array_2D(self, array, new_size):
        """Crop array to smaller size by taking center portion"""
        old_size = array.shape

        # Calculate offsets
        offset_x = (old_size[0] - new_size) // 2
        offset_y = (old_size[1] - new_size) // 2

        # Extract center portion
        return array[offset_x:offset_x + new_size, offset_y:offset_y + new_size]


class dipdistr:
    """
    Calculates the theoretical point spread function (PSF) for fixed dipoles.

    This version uses the approach from combined_hinterer.m converted to Python,
    while maintaining the same interface as the original dipdistr class.
    """

    def __init__(self, wl, n, n0, M, NA):
        # Original parameters
        self.wl = wl  # Wavelength in nm (emission peak in buffer)
        self.NA = NA  # Numerical aperture
        self.M = M  # Magnification
        self.n = n  # Refractive index of glass/oil
        self.np = 1.0  # Refractive index of air in lab (vacuum)
        self.n0 = n0  # Refractive index of sample medium (water)

        self.kp = 2 * np.pi / wl  # Wavevector amplitude in air (vacuum)
        self.k0 = n0 * self.kp  # Wavevector amplitude in buffer
        self.k = n * self.kp  # Wavevector amplitude in glass

        self.etapmed = n0 / M  # Integration limits
        self.etapmax = NA / M

        # Current position and focus
        self.rho = 0  # Current radial position
        self.deltaz = 0  # Current defocus position

        # Calculate or load normalization constants
        try:
            self.norm = np.loadtxt('dipoletablenorm.dat')
        except:
            self.deltaz = -1.0
            self.norm = self.Normalization()
            np.savetxt('dipoletablenorm.dat', self.norm)

        # Keep track of parameters for tabulation
        self.focusvals = np.arange(-150, 150, 1.0)
        self.nfoci = len(self.focusvals)
        self.focusindex = 0  # Will be set when calculating PSF

        # Parameters from MATLAB implementation
        self.n_discretization_bfp = 129
        self.n_pixels = 170
        self.oversampling = 9
        self.pixel_size = 51.2e-10  # In meters
        self.unit_object_space = self.pixel_size / self.oversampling
        self.unit_k_space = 2 * self.NA / (self.wl * 1e-9) * 2 * np.pi / self.n_discretization_bfp
        self.objective_focal_length = 3e-3  # In meters

        # Cache for expensive computations
        self.psf_cache = {}
        self.bfp_cache = {}

    # --- Methods to maintain compatibility with original API ---

    def SetFocus(self, deltaz):
        """Set the focus (defocus distance)"""
        self.deltaz = deltaz
        # Find closest tabulated focus value for compatibility
        self.focusindex = np.argmin(np.abs(self.focusvals - deltaz))

    def GetFocus(self):
        """Get the current focus setting"""
        return self.deltaz

    def SetRho(self, rho):
        """Set the radial position"""
        self.rho = rho

    def GetRho(self):
        """Get the current radial position"""
        return self.rho

    def Normalization(self):
        """Calculate normalization constants for different orientations"""
        print("Calculating normalization constants...")
        # Simplified estimation: calculate for in-plane and out-of-plane dipoles
        return (1.0, 1.0)

    def _create_psf_object(self, x, y, alpha, beta, deltaz):
        """Create a PSF object with the given parameters"""

        class PsfObject:
            """Simple container class to hold PSF parameters"""
            pass

        psf = PsfObject()

        # Physical parameters
        psf.wavelength = self.wl * 1e-9  # Convert to meters
        psf.refractive_indices = [self.np, self.n, self.n0]  # [specimen, intermediate, immersion]

        # Microscope parameters
        psf.n_discretization_bfp = self.n_discretization_bfp
        psf.n_pixels = self.n_pixels
        psf.oversampling = self.oversampling
        psf.unit_object_space = self.unit_object_space
        psf.unit_k_space = self.unit_k_space
        psf.objective_focal_length = self.objective_focal_length

        # Convert coordinates to meters
        psf.position = np.array([x, y, deltaz]) * 1e-9
        psf.height_intermediate_layer = 0

        # Setup dipole with given angles
        psf.dipole = Dipole(beta, alpha)

        # Setup empty masks
        psf.pupil_mask = Mask(psf.n_discretization_bfp)

        return psf

    def _calculate_aberrations(self, psf_obj):
        """Calculate aberrations for the given PSF object"""
        # Get Zernike polynomials calculator
        zernike = ZernikePolynomials.get_instance(psf_obj.pupil_mask)

        # Calculate tilt aberrations from position
        x, y = psf_obj.position[0:2]
        unit_k_space = psf_obj.unit_k_space
        unit_object_space = psf_obj.unit_object_space
        n_bfp = psf_obj.n_discretization_bfp

        # Convert position to coefficients for tilt Zernike polynomials (indices 2 and 3)
        coeff_factor = unit_object_space * unit_k_space * n_bfp / (4 * np.pi)
        coeff_x = -x * coeff_factor
        coeff_y = -y * coeff_factor

        # Calculate aberrations (tilt only for now)
        # Note: We keep using indices 2 and 3 as they correspond to x and y tilts in Noll indexing
        aberrations = zernike.get_aberration([2, 3], [coeff_x, coeff_y])

        # Add defocus (Zernike index 4 corresponds to defocus)
        if psf_obj.position[2] != 0:
            # Calculate defocus coefficient
            z = psf_obj.position[2]
            wavelength = psf_obj.wavelength
            # Simplified formula for defocus coefficient
            defocus_coeff = z / (4 * wavelength)

            # Add defocus term (index 4)
            defocus_aberration = zernike.get_aberration([4], [defocus_coeff])
            aberrations += defocus_aberration

        return aberrations

    def _apply_aberrations(self, psf_obj, E_BFP, aberrations):
        """Apply aberrations to the electric field"""
        # Apply phase shift from aberrations
        mask = psf_obj.pupil_mask.values * np.exp(1j * 2 * np.pi * aberrations)

        # Apply to both field components
        E_BFP_aberrated = {
            'x': E_BFP['x'] * mask,
            'y': E_BFP['y'] * mask
        }

        return E_BFP_aberrated

    def _calculate_psf(self, x, y, alpha, beta, deltaz):
        """Calculate the PSF using the methods from MATLAB"""
        # Create a PSF object
        psf_obj = self._create_psf_object(x, y, alpha, beta, deltaz)

        # Calculate electric field in back focal plane
        bfp = BackFocalPlane(psf_obj)
        E_BFP = bfp.electric_field

        # Calculate aberrations
        aberrations = self._calculate_aberrations(psf_obj)

        # Apply aberrations
        E_BFP_aberrated = self._apply_aberrations(psf_obj, E_BFP, aberrations)

        # Initialize ChirpZ transform
        chirp_z = ChirpZTransform(psf_obj)

        # Transform to image plane
        E_image_x = chirp_z.apply(psf_obj, E_BFP_aberrated['x'])
        E_image_y = chirp_z.apply(psf_obj, E_BFP_aberrated['y'])

        # Calculate intensities
        Ix = np.abs(E_image_x) ** 2
        Iy = np.abs(E_image_y) ** 2

        # Total intensity
        I_total = Ix + Iy

        # Normalization
        total_intensity = np.sum(I_total)
        if total_intensity > 0:
            I_total = I_total / total_intensity

        # fig, ax = plt.subplots(figsize=(10, 10))  # No need for (1,1), just use `ax`
        # im = ax.imshow(np.abs(E_BFP['y'])**2, cmap='gray')
        # ax.set_xlabel('x (nm)')
        # ax.set_ylabel('y (nm)')
        # ax.grid(False)
        # plt.tight_layout()
        # plt.show()
        # plt.close()

        return I_total

    def _get_intensity_at_point(self, psf_image, r_nm, phi):
        """Extract intensity at a specific point in the PSF image"""
        # Convert r_nm and phi to pixel coordinates
        center = psf_image.shape[0] // 2

        # Fix: Handle even/odd array sizes correctly
        offset = 0 if psf_image.shape[0] % 2 == 1 else 0.5

        # Convert to pixels
        r_pixels = r_nm / (self.unit_object_space * 1e9)

        # Convert to Cartesian (using correct arctan2)
        x_pixels = int(round(r_pixels * np.cos(phi) + offset))
        y_pixels = int(round(r_pixels * np.sin(phi) + offset))

        # Get intensity at that pixel (if within range)
        if (0 <= center + y_pixels < psf_image.shape[0] and
                0 <= center + x_pixels < psf_image.shape[1]):
            return psf_image[center + y_pixels, center + x_pixels]
        else:
            # Outside image bounds, return very low intensity
            return 1e-10

    # --- Main PSF calculation methods ---

    def PSF_exact(self, x, y, alpha, beta, deltaz):
        """
        Calculate the exact PSF using the BackFocalPlane approach.
        This maintains the same interface as the original method.

        Parameters:
        x, y: coordinates in the image plane (nm)
        alpha: azimuthal angle of the dipole (rad)
        beta: polar angle of the dipole (rad)
        deltaz: defocus distance (nm)

        Returns:
        Intensity value at the specified point
        """
        # Fix: Make alpha transformation consistent with Zernike conventions
        # We keep this transformation for backward compatibility
        alpha_adjusted = alpha + np.pi

        # Set current defocus
        self.SetFocus(deltaz)

        # # Cache key for this orientation and focus
        # cache_key = f"{beta:.4f}_{alpha_adjusted:.4f}_{deltaz:.1f}"
        #
        # # Calculate or retrieve PSF image
        # if cache_key not in self.psf_cache:
        #     self.psf_cache[cache_key] = self._calculate_psf(0, 0, alpha_adjusted, beta, deltaz)
        #
        # # Get the full PSF image
        # psf_image = self.psf_cache[cache_key]

        # No caching version
        psf_image = self._calculate_psf(0, 0, alpha_adjusted, beta, deltaz)


        # Convert (x,y) to polar coordinates using arctan2
        r = np.sqrt(x ** 2 + y ** 2)
        phi = np.arctan2(y, x)

        # Get the intensity value at the specific point
        intensity = self._get_intensity_at_point(psf_image, r, phi)

        # Scale by magnification squared as in original
        return self.M ** 2 * intensity

    def PSF_approx(self, x, y, alpha, beta, deltaz):
        """
        Calculate the approximate PSF.
        For consistency with the original API, this calls the same method as PSF_exact.
        """
        return self.PSF_exact(x, y, alpha, beta, deltaz)