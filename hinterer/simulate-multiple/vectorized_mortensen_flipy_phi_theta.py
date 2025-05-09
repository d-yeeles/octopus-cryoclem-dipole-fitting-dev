import numpy as np
from MLEwT_fixed import dipdistr
from math import ceil
import os
import time
from scipy.special import jn

class VectorizedDipolePSFGenerator:
    def __init__(self, image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor=3, verbose=False):
        self.image_size = (int(image_size[0]), int(image_size[1]))
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
        
        # Initialize dipole distribution model
        self.DD = dipdistr(wavelength, n_objective, n_sample, magnification, NA, norm_file)
    
    def __call__(self, phi, theta, x_pos, y_pos, n_photons, verbose=False):
        """
        Vectorized implementation of the PSF generation
        """
        # To match Hinterer
        phi = -phi
        theta = -theta

        # Create position vectors
        y_pos_vec = np.arange(-(self.oversampled_image_size[0]-1)/2, self.oversampled_image_size[0]/2) * self.oversampled_pixel_size
        x_pos_vec = np.arange(-(self.oversampled_image_size[1]-1)/2, self.oversampled_image_size[1]/2) * self.oversampled_pixel_size
        
        # Create 2D meshgrids for the coordinates
        X, Y = np.meshgrid(x_pos_vec, y_pos_vec)
        
        # Adjust for the dipole position
        X = X - x_pos
        Y = -Y - y_pos # To match Hinterer
        
        # Calculate radial distance and angle
        r = np.sqrt(X**2 + Y**2)
        phip = np.arctan2(Y, X) % (2*np.pi)
        
        # Calculate scaled radius (rho)
        rho = self.DD.M * r
        
        # Initialize PSF array
        dipole_psf = np.zeros(self.oversampled_image_size)
        
        # Use the existing DD.PSF_approx but in a more efficient way
        for i in range(self.oversampled_image_size[0]):
            for j in range(self.oversampled_image_size[1]):
                dipole_psf[i, j] = self.DD.PSF_approx(X[i, j], Y[i, j], phi, theta)

        # Reshape for subpixel averaging
        dipole_psf = dipole_psf.reshape(self.image_size[0], self.subpixel_factor, 
                                     self.image_size[1], self.subpixel_factor)
        dipole_psf = np.mean(dipole_psf, axis=(1, 3))



        # Mask off central 80nm radius circle patch to avoid glowing corners
        # Create and apply mask with a specific background value
        radius_nm = 500  # Set fixed 520nm radius
        background_value = 0.0000001  # Replace 0.1 with your desired background value

        # Create coordinates for the final image
        y_final = np.arange(-(self.image_size[0]-1)/2, self.image_size[0]/2) * self.pixel_size
        x_final = np.arange(-(self.image_size[1]-1)/2, self.image_size[1]/2) * self.pixel_size
        X_final, Y_final = np.meshgrid(x_final, y_final)

        # Calculate distances from dipole location
        X_final = X_final - x_pos
        Y_final = -Y_final - y_pos  # To match Hinterer
        r_final = np.sqrt(X_final**2 + Y_final**2)

        # Create mask (1 inside radius, 0 outside)
        mask = (r_final <= radius_nm)

        # Apply mask: keep original values inside circle, set to background_value outside
        dipole_psf = np.where(mask, dipole_psf, background_value)




        
        # Scale to the correct number of photons
        dipole_psf = dipole_psf / np.sum(dipole_psf) * n_photons
        
        return dipole_psf


class FullyVectorizedDipolePSFGenerator:
    """
    A more comprehensive vectorization approach that recreates the PSF_approx functionality
    in a fully vectorized manner, eliminating the need for pixel-by-pixel processing.
    """
    def __init__(self, image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file, subpixel_factor=3, verbose=False):
        self.image_size = (int(image_size[0]), int(image_size[1]))
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
        
        # Initialize dipole distribution model (we'll use this for parameters only)
        self.DD = dipdistr(wavelength, n_objective, n_sample, magnification, NA, norm_file)
    
    def _vectorized_bessel0cum_subcrit(self, rho):
        """Vectorized implementation of Bessel0Cum_subcrit"""
        kp = self.DD.kp
        norm = self.DD.b0sub_norm
        mom1 = self.DD.b0sub_mean
        mom2 = self.DD.b0sub_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j0 = jn(0, arg[mask_nonsmall])
            j1 = jn(1, arg[mask_nonsmall])
            j2 = jn(2, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j0 + 0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * (j2 - j1/arg[mask_nonsmall]))
        
        return result
    
    def _vectorized_bessel0cum_real(self, rho):
        """Vectorized implementation of Bessel0Cum_real"""
        kp = self.DD.kp
        norm = self.DD.b0real_norm
        mom1 = self.DD.b0real_mean
        mom2 = self.DD.b0real_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            result[mask_nonsmall] = norm * (jn(0, arg[mask_nonsmall]) + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (jn(2, arg[mask_nonsmall]) - jn(1, arg[mask_nonsmall])/arg[mask_nonsmall]))
        
        return result
    
    def _vectorized_bessel0cum_imag(self, rho):
        """Vectorized implementation of Bessel0Cum_imag"""
        kp = self.DD.kp
        norm = self.DD.b0imag_norm
        mom1 = self.DD.b0imag_mean
        mom2 = self.DD.b0imag_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            result[mask_nonsmall] = norm * (jn(0, arg[mask_nonsmall]) + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (jn(2, arg[mask_nonsmall]) - jn(1, arg[mask_nonsmall])/arg[mask_nonsmall]))
        
        return result
    
    def _vectorized_bessel1cum_subcrit(self, rho):
        """Vectorized implementation of Bessel1Cum_subcrit"""
        kp = self.DD.kp
        norm = self.DD.b1sub_norm
        mom1 = self.DD.b1sub_mean
        mom2 = self.DD.b1sub_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j1 = jn(1, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j1 + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (2 * j1 / arg[mask_nonsmall]**2 - 
                                            jn(0, arg[mask_nonsmall]) / arg[mask_nonsmall] - j1))
        
        return result
    
    def _vectorized_bessel1cum_real(self, rho):
        """Vectorized implementation of Bessel1Cum_real"""
        kp = self.DD.kp
        norm = self.DD.b1real_norm
        mom1 = self.DD.b1real_mean
        mom2 = self.DD.b1real_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j1 = jn(1, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j1 + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (2 * j1 / arg[mask_nonsmall]**2 - 
                                            jn(0, arg[mask_nonsmall]) / arg[mask_nonsmall] - j1))
        
        return result
    
    def _vectorized_bessel1cum_imag(self, rho):
        """Vectorized implementation of Bessel1Cum_imag"""
        kp = self.DD.kp
        norm = self.DD.b1imag_norm
        mom1 = self.DD.b1imag_mean
        mom2 = self.DD.b1imag_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j1 = jn(1, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j1 + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (2 * j1 / arg[mask_nonsmall]**2 - 
                                            jn(0, arg[mask_nonsmall]) / arg[mask_nonsmall] - j1))
        
        return result
    
    def _vectorized_bessel2cum_subcrit(self, rho):
        """Vectorized implementation of Bessel2Cum_subcrit"""
        kp = self.DD.kp
        norm = self.DD.b2sub_norm
        mom1 = self.DD.b2sub_mean
        mom2 = self.DD.b2sub_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j2 = jn(2, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j2 + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (jn(0, arg[mask_nonsmall]) - 
                                            3 * jn(1, arg[mask_nonsmall]) / arg[mask_nonsmall] + 
                                            6 * j2 / arg[mask_nonsmall]**2))
        
        return result
    
    def _vectorized_bessel2cum_real(self, rho):
        """Vectorized implementation of Bessel2Cum_real"""
        kp = self.DD.kp
        norm = self.DD.b2real_norm
        mom1 = self.DD.b2real_mean
        mom2 = self.DD.b2real_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j2 = jn(2, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j2 + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (jn(0, arg[mask_nonsmall]) - 
                                            3 * jn(1, arg[mask_nonsmall]) / arg[mask_nonsmall] + 
                                            6 * j2 / arg[mask_nonsmall]**2))
        
        return result
    
    def _vectorized_bessel2cum_imag(self, rho):
        """Vectorized implementation of Bessel2Cum_imag"""
        kp = self.DD.kp
        norm = self.DD.b2imag_norm
        mom1 = self.DD.b2imag_mean
        mom2 = self.DD.b2imag_var
        
        arg = kp * rho * mom1
        
        # Handle case where arg is zero or very small
        mask_small = np.abs(arg) < 1e-10
        result = np.zeros_like(rho)
        
        # For small args, just return norm
        result[mask_small] = norm
        
        # For non-small args, calculate the value
        mask_nonsmall = ~mask_small
        if np.any(mask_nonsmall):
            j2 = jn(2, arg[mask_nonsmall])
            result[mask_nonsmall] = norm * (j2 + 
                                           0.5 * (rho[mask_nonsmall] * kp)**2 * mom2 * 
                                           (jn(0, arg[mask_nonsmall]) - 
                                            3 * jn(1, arg[mask_nonsmall]) / arg[mask_nonsmall] + 
                                            6 * j2 / arg[mask_nonsmall]**2))
        
        return result
    
    def _vectorized_intensity_approx(self, rho, phip, alpha, beta):

        """Vectorized implementation of Intensity_approx"""
        # Calculate all the Bessel function values
        bessel0cum_subcrit = self._vectorized_bessel0cum_subcrit(rho)
        bessel0cum_real = self._vectorized_bessel0cum_real(rho)
        bessel0cum_imag = self._vectorized_bessel0cum_imag(rho)
        
        bessel1cum_subcrit = self._vectorized_bessel1cum_subcrit(rho)
        bessel1cum_real = self._vectorized_bessel1cum_real(rho)
        bessel1cum_imag = self._vectorized_bessel1cum_imag(rho)
        
        bessel2cum_subcrit = self._vectorized_bessel2cum_subcrit(rho)
        bessel2cum_real = self._vectorized_bessel2cum_real(rho)
        bessel2cum_imag = self._vectorized_bessel2cum_imag(rho)
        
        # Calculate the final intensity
        sin_beta = np.sin(beta)
        cos_beta = np.cos(beta)
        
        par = sin_beta**2 / 4.0 * (
            (bessel0cum_subcrit + bessel0cum_real)**2 + bessel0cum_imag**2 +
            (bessel2cum_subcrit + bessel2cum_real)**2 + bessel2cum_imag**2 -
            2 * np.cos(2 * (phip - alpha)) * (
                (bessel0cum_subcrit + bessel0cum_real) *
                (bessel2cum_subcrit + bessel2cum_real) +
                bessel0cum_imag * bessel2cum_imag
            )
        )
        
        mix = sin_beta * cos_beta * np.cos(phip - alpha) * (
            (bessel1cum_subcrit + bessel1cum_real) * bessel2cum_imag -
            (bessel2cum_subcrit + bessel2cum_real) * bessel1cum_imag -
            (bessel1cum_subcrit + bessel1cum_real) * bessel0cum_imag +
            (bessel0cum_subcrit + bessel0cum_real) * bessel1cum_imag
        )
        
        vert = cos_beta**2 * (
            (bessel1cum_subcrit + bessel1cum_real)**2 + bessel1cum_imag**2
        )
        
        value = par + mix + vert
        
        # Normalization
        value /= (sin_beta**2 * self.DD.norm[0] + cos_beta**2 * self.DD.norm[1])
        
        return value
    
    def __call__(self, phi, theta, x_pos, y_pos, n_photons, verbose=False):
        """
        Fully vectorized PSF generation
        """
        start_time = time.time()
       
        # To match Hinterer
        phi = -phi
        theta = -theta
 
        # Create meshgrid for position vectors
        y_pos_vec = np.arange(-(self.oversampled_image_size[0]-1)/2, self.oversampled_image_size[0]/2) * self.oversampled_pixel_size
        x_pos_vec = np.arange(-(self.oversampled_image_size[1]-1)/2, self.oversampled_image_size[1]/2) * self.oversampled_pixel_size
        
        # Create 2D meshgrids for the coordinates
        X, Y = np.meshgrid(x_pos_vec, y_pos_vec)
        
        # Adjust for the dipole position
        X = X - x_pos
        Y = -Y - y_pos # To match Hinterer
        
        # Calculate radial distance and angle
        r = np.sqrt(X**2 + Y**2)
        phip = np.arctan2(Y, X) % (2*np.pi)
        
        # Calculate scaled radius (rho)
        rho = self.DD.M * r
        
        # Calculate PSF using vectorized intensity method
        intensity = self._vectorized_intensity_approx(rho, phip, phi, theta)
        dipole_psf = self.DD.M**2 * intensity
        
        # Reshape for subpixel averaging
        dipole_psf = dipole_psf.reshape(self.image_size[0], self.subpixel_factor, 
                                      self.image_size[1], self.subpixel_factor)
        dipole_psf = np.mean(dipole_psf, axis=(1, 3))
       




        # Mask off central 80nm radius circle patch to avoid glowing corners
        # Create and apply mask with a specific background value
        radius_nm = 500  # Set fixed 520nm radius
        background_value = 0.0000001  # Replace 0.1 with your desired background value

        # Create coordinates for the final image
        y_final = np.arange(-(self.image_size[0]-1)/2, self.image_size[0]/2) * self.pixel_size
        x_final = np.arange(-(self.image_size[1]-1)/2, self.image_size[1]/2) * self.pixel_size
        X_final, Y_final = np.meshgrid(x_final, y_final)

        # Calculate distances from dipole location
        X_final = X_final - x_pos
        Y_final = -Y_final - y_pos  # To match Hinterer
        r_final = np.sqrt(X_final**2 + Y_final**2)

        # Create mask (1 inside radius, 0 outside)
        mask = (r_final <= radius_nm)

        # Apply mask: keep original values inside circle, set to background_value outside
        dipole_psf = np.where(mask, dipole_psf, background_value)





 
        # Scale to the correct number of photons
        dipole_psf = dipole_psf / np.sum(dipole_psf) * n_photons
        
        end_time = time.time()
#        if verbose:
#            print(f"Vectorized PSF generation took {end_time - start_time:.4f} seconds")
        
        return dipole_psf


def run_simulator_vectorized(x, y, theta, phi, image_size_px, pixel_size, wavelength, n_objective, n_sample, NA, n_photons):
    """
    Runs the vectorized Mortensen simulator.
    """
    # Define parameters
    image_size = (image_size_px, image_size_px)  
    magnification = 215 
    norm_file = os.path.join(os.path.expanduser('~'), 'dipolenorm.npy')
    
    # Factor by which to oversample
    subpixel_factor = 1
    
    # Create PSF generator instance
    psf_generator = FullyVectorizedDipolePSFGenerator(
        image_size, pixel_size, wavelength, n_objective, n_sample, 
        magnification, NA, norm_file, subpixel_factor, verbose=True
    )
    
    # Generate dipole PSF
    dipole_psf = psf_generator(phi, theta, x, y, n_photons, verbose=True)
    
    return dipole_psf


def benchmark_comparison():
    """
    Compare the original and vectorized implementations.
    """
    # Define test parameters
    image_size_px = 21
    pixel_size = 80
    wavelength = 680
    n_objective = 1.518
    n_sample = 1.33
    NA = 1.49
    n_photons = 1000
    x = 0
    y = 0
    theta = 0.5
    phi = 1.0
    
    try:
        # Original implementation
        from mortensen_simulator import DipolePSFGenerator, run_simulator
        
        print("Testing original implementation...")
        start_time = time.time()
        original_psf = run_simulator(x, y, theta, phi, image_size_px, pixel_size, wavelength, n_objective, n_sample, NA, n_photons)
        original_time = time.time() - start_time
        print(f"Original implementation took: {original_time:.4f} seconds")
        
        # Vectorized implementation (partially vectorized)
        print("\nTesting partially vectorized implementation...")
        psf_generator1 = VectorizedDipolePSFGenerator(
            (image_size_px, image_size_px), pixel_size, wavelength, n_objective, n_sample, 
            215, NA, os.path.join(os.path.expanduser('~'), 'dipolenorm.npy'), 9
        )
        start_time = time.time()
        partial_psf = psf_generator1(phi, theta, x, y, n_photons)
        partial_time = time.time() - start_time
        print(f"Partially vectorized implementation took: {partial_time:.4f} seconds")
        print(f"Speedup vs original: {original_time/partial_time:.2f}x")
        
        # Fully vectorized implementation
        print("\nTesting fully vectorized implementation...")
        start_time = time.time()
        vectorized_psf = run_simulator_vectorized(x, y, theta, phi, image_size_px, pixel_size, wavelength, n_objective, n_sample, NA, n_photons)
        vectorized_time = time.time() - start_time
        print(f"Fully vectorized implementation took: {vectorized_time:.4f} seconds")
        print(f"Speedup vs original: {original_time/vectorized_time:.2f}x")
        
        # Compare results
        diff1 = np.abs(original_psf - partial_psf)
        diff2 = np.abs(original_psf - vectorized_psf)
        print("\nComparison to original implementation:")
        print(f"Partially vectorized - Max difference: {np.max(diff1):.6f}, Mean difference: {np.mean(diff1):.6f}")
        print(f"Fully vectorized - Max difference: {np.max(diff2):.6f}, Mean difference: {np.mean(diff2):.6f}")
        
        return original_psf, partial_psf, vectorized_psf
        
    except ImportError:
        print("Original mortensen_simulator.py module not found. Skipping comparison.")
        
        # Just run the vectorized implementation
        start_time = time.time()
        vectorized_psf = run_simulator_vectorized(x, y, theta, phi, image_size_px, pixel_size, wavelength, n_objective, n_sample, NA, n_photons)
        vectorized_time = time.time() - start_time
        print(f"Fully vectorized implementation took: {vectorized_time:.4f} seconds")
        
        return None, None, vectorized_psf


if __name__ == "__main__":
    print("Testing vectorized Mortensen PSF simulator...")
    results = benchmark_comparison()
