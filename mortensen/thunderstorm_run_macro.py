import imagej
import time
import pandas as pd

import os
os.environ['JAVA_HOME'] = '/usr/lib/jvm/java-11-openjdk-amd64/'

# Initialize ImageJ (this will download Fiji if not already installed)
ij = imagej.init('/home/tfq96423/fiji-linux64/Fiji.app', headless=False)

def run_thunderstorm(frames_path, results_path):

    # Function to run ImageJ macro
    def run_imagej_macro(macro_code):
        # Run the macro
        ij.ui().showUI()  # Ensures graphical UI is displayed if needed
        ij.py.run_macro(macro_code)

        # Add time to allow graphical output to render
        time.sleep(3)  # Adjust this as needed

    # Load your macro (adjust the path if needed)
    macro_code = f"""
      File.openSequence("{frames_path}");
      run("Camera setup", "offset=414.0 isemgain=false photons2adu=3.6 pixelsize=23.4");
      run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.6 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false renderer=[Averaged shifted histograms] magnification=5.0 colorizez=false threed=false shifts=2 repaint=50");
      run("Export results", "filepath={results_path} fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=false x=true y=true bkgstd=true id=false uncertainty=true frame=true");
      close();
    """

    # Call the function
    run_imagej_macro(macro_code)

    # Close ImageJ properly after execution
    ij.dispose()

    # # Import the positions in each frame
    # df = pd.read_csv(results_path)
    # df_selected_columns = df[['frame','x [nm]', 'y [nm]']]
    # f_array = df['frame'].to_numpy()  # 1D array for frame numbers
    # x_array = df['x [nm]'].to_numpy()  # 1D array for x values
    # y_array = df['y [nm]'].to_numpy()  # 1D array for y values

    return# f_array, x_array, y_array

def reconstruct(results_path, output_img_path):

    # Function to run ImageJ macro
    def run_imagej_macro(macro_code):
        # Run the macro
        ij.ui().showUI()  # Ensures graphical UI is displayed if needed
        ij.py.run_macro(macro_code)

        # Add time to allow graphical output to render
        time.sleep(3)  # Adjust this as needed


    macro_code = f"""
        run("Import results", "filepath={results_path} fileformat=[CSV (comma separated)] livepreview=true rawimagestack= startingframe=1 append=false");
        selectImage("Averaged shifted histograms");
        run("Image...  ", "outputfile={output_img_path}");
    """

    # Call the function
    run_imagej_macro(macro_code)

    # Close ImageJ properly after execution
    ij.dispose()

    return