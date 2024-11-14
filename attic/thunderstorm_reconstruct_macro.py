import imagej
import time


import os
os.environ['JAVA_HOME'] = '/usr/lib/jvm/java-11-openjdk-amd64/'

# Initialize ImageJ (this will download Fiji if not already installed)
ij_new = imagej.init('/home/tfq96423/fiji-linux64/Fiji.app', headless=False)

def reconstruct(results_path, output_img_path):

    # Function to run ImageJ macro
    def run_imagej_macro(macro_code):
        # Run the macro
        ij_new.ui().showUI()  # Ensures graphical UI is displayed if needed
        ij_new.py.run_macro(macro_code)

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
    ij_new.dispose()

    return