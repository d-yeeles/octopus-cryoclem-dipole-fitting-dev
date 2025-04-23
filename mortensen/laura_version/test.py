import numpy as np
import sys
import datetime
import matplotlib.pyplot as plt
#from MLEwT_fixed import dipdistr, MLEwT
from MLEwT_matlab_engine import dipdistr, MLEwT

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

    def mortensen_fit(self, dipole_psf, init_theta, init_phi):
        
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

def run_mortensen_fit(phi, theta):
    """
    Runs the Mortensen fit for given phi and theta, returning the results and ground truth.
    """
    # Define parameters
    image_size = (18, 18)  
    pixel_size = 51  
    wavelength = 500  
    n_objective = 2.17  
    n_sample = 1.31  
    magnification = 215  
    NA = 2.17  
    norm_file = "/home/tfq96423/dipolenorm.npy"
    n_photons = 2000

    # Define dipole ground_truth position in nm
    x_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2)  
    y_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2) 
    
    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file)    

    # Generate dipole PSF
    dipole_psf = psf_generator(phi, theta, x_pos_nm, y_pos_nm, n_photons)
    dipole_psf_noisy = np.random.poisson(dipole_psf)

    # Run Mortensen fit
    results = psf_generator.mortensen_fit(dipole_psf_noisy, theta, phi)
    
    # Return results instead of printing
    return results, [phi, theta, x_pos_nm, y_pos_nm, n_photons]

# Prevent script from running when imported
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python test.py <phi> <theta>")
        sys.exit(1)

    phi = float(sys.argv[1])
    theta = float(sys.argv[2])

    results, ground_truth = run_mortensen_fit(phi, theta)
    
#    print(f"Results from the Mortensen fit are:  {', '.join(map(str, results))}")    
#    print(f"Ground truth are: {', '.join(map(str, ground_truth))}")
