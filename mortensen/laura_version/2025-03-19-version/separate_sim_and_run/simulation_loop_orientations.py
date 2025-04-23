import numpy as np
from MLEwT_fixed import dipdistr, MLEwT
import tifffile as tiff

def save_as_tif(psf_total_image, output_path):
    # Ensure the image is uint32
    psf_total_image = psf_total_image.astype(np.uint32)
    
    # Save as TIF with LZW compression
    with tiff.TiffWriter(output_path, bigtiff=True) as t:
        t.write(
            psf_total_image, 
            photometric='minisblack',  # Grayscale
            dtype=np.uint32,  # 32-bit integer
            rowsperstrip=16,  # Strip length
            compression='lzw',  # Lossless compression
            software='Python'  # Software metadata
        )


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


def run_simulator(phi, theta):
    """
    Runs the Mortensen fit for given phi and theta, returning the results and ground truth.
    """
    # Define parameters
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

    # Define dipole ground_truth position in nm
    x_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2)  
    y_pos_nm = np.random.uniform(0 - pixel_size/2, 0 + pixel_size/2) 
    
    # Create PSF generator instance
    psf_generator = DipolePSFGenerator(image_size, pixel_size, wavelength, n_objective, n_sample, magnification, NA, norm_file)    

    # Generate dipole PSF
    dipole_psf = psf_generator(phi, theta, x_pos_nm, y_pos_nm, n_photons)
    dipole_psf_noisy = np.random.poisson(dipole_psf)

    parameters = [number_of_spots, pixel_size, pixel_size*image_size[0], image_size[0], wavelength, NA, objectiveFocalLength, [n_sample, n_objective, n_objective], nDiscretizationBFP, background, n_photons, x_pos_nm, y_pos_nm, theta, phi]

    # Return results instead of printing
    return dipole_psf_noisy, parameters

# Prevent script from running when imported
if __name__ == "__main__":

    runs = range(1)
    thetas = [67.5*np.pi/180]#np.arange(0, np.pi/2 + 0.001, 22.5*np.pi/180)
    phis = np.arange(0, 2*np.pi, 45*np.pi/180)

    for run in runs:
      for theta in thetas:
        for phi in phis:

          theta_deg = theta*180/np.pi
          phi_deg = phi*180/np.pi

          psf_image, parameters = run_simulator(phi, theta)
          save_as_tif(psf_image, f'./sims/sim_theta{round(theta_deg):03}_phi{round(phi_deg):03}_run{int(run+1)}.tif')


          # Define output path
          params_output_path = f"./sims/params_inc{round(theta_deg):03}_az{round(phi_deg):03}_run{int(run+1)}.m"
    
          # Write data to file
          with open(params_output_path, 'w') as file:
            file.write(f"% ground truth for sim_inc{round(theta_deg)}_az{round(phi_deg)}_run{int(run+1)}.tif\n")
            file.write("% settings\n")
            file.write(f"number_of_spots = {parameters[0]}\n")
            file.write(f"pixel_size_nm = {parameters[1]}\n")
            file.write(f"image_size_nm = {parameters[2]}\n")
            file.write(f"image_size_px = {parameters[3]}\n")
            file.write(f"wavelength = {parameters[4]}\n")
            file.write(f"par.objectiveNA = {parameters[5]}\n")
            file.write(f"objectiveFocalLength = {parameters[6]}\n")
            file.write(f"par.refractiveIndices = {parameters[7]}\n")
            file.write(f"par.nDiscretizationBFP = {parameters[8]}\n")
            file.write(f"par.backgroundNoise = {parameters[9]}\n")
            file.write(f"par.nPhotons = {parameters[10]}\n")
            file.write(f"positionX_nm_array = [{parameters[11]}, ]\n")
            file.write(f"positionY_nm_array = [{parameters[12]}, ]\n")
            file.write(f"angleInclination_array = [{parameters[13]}, ]\n")
            file.write(f"angleAzimuth_array = [{parameters[14]}, ]\n")


