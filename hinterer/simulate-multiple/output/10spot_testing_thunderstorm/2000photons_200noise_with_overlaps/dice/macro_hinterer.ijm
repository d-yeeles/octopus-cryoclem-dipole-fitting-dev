// Base path for export
basePath = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/10spot_testing_thunderstorm/2000photons_200noise_with_overlaps/dice/results/";

// Create a log file to track progress
logFile = basePath + "processing_log.txt";
File.saveString("Starting parameter sweep at " + getDateTime() + "\n", logFile);

// Loop over sigma values from 1 to 4 in steps of 0.25
for (sigma = 1.0; sigma <= 4.0; sigma += 0.25) {
    // Format sigma to 2 decimal places to avoid floating point precision issues
    sigmaStr = d2s(sigma, 2);
    
    // Loop over fitradius values from 1 to 6 in steps of 1
    for (fitradius = 1; fitradius <= 6; fitradius++) {
        // Log current parameters
        logMessage = "Processing: sigma = " + sigmaStr + ", fitradius = " + fitradius + " at " + getDateTime() + "\n";
        File.append(logMessage, logFile);

        // Add a longer delay to allow GUI to update properly
        wait(1000);
        
        // Run the analysis with current parameters
        run("Run analysis", "filter=[Lowered Gaussian filter] sigma=" + sigmaStr + 
            " detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) " + 
            "estimator=[PSF: Gaussian] sigma=1.6 fitradius=" + fitradius + 
            " method=[Maximum likelihood] full_image_fitting=false mfaenabled=false renderer=[No Renderer]");

        // Add a longer delay to allow GUI to update properly
        wait(2000);
        
        
        // Create filename reflecting the current parameter values
        filename = "thunderstorm-hinterer-sigma" + sigmaStr + "-fitrad" + fitradius + ".csv";
        
        // Export results with the appropriate filename
        run("Export results", "filepath=" + basePath + filename + 
            " fileformat=[CSV (comma separated)] sigma=true intensity=true offset=true " + 
            "saveprotocol=false x=true y=true bkgstd=true id=true uncertainty=true frame=true");
        
        // Close the results table to prevent GUI issues
        if (isOpen("Results")) {
            selectWindow("Results");
            run("Close");
        }
        
        // Add a longer delay to allow GUI to update properly
        wait(1000);
    
        // Force garbage collection between major iterations
        call("java.lang.System.gc");

    }
}

// Log completion
File.append("Processing complete at " + getDateTime() + "\n", logFile);

// Function to get current date and time
function getDateTime() {
    getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
    return "" + year + "-" + IJ.pad(month+1, 2) + "-" + IJ.pad(dayOfMonth, 2) + " " + 
           IJ.pad(hour, 2) + ":" + IJ.pad(minute, 2) + ":" + IJ.pad(second, 2);
}
