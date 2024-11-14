run("Import results", "filepath=/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/mortensen_results.csv fileformat=[CSV (comma separated)] livepreview=true rawimagestack= startingframe=1 append=false");

wait(3000);

saveAs("Tiff", "/home/tfq96423/Documents/cryoCLEM/dipole-issue/mortensen-loop/reconstruction.tif");

close();
