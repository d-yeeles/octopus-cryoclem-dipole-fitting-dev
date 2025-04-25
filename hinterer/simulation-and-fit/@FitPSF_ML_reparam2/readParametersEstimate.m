function par = readParametersEstimate(psf)

    % dave apr 2025 - add this so we can inherit the right params when
    % fitting
    par.nPixels = psf.nPixels;
    par.pixelSize = psf.pixelSize;
    par.wavelength = psf.wavelength;
    par.objectiveNA = psf.objectiveNA;
    par.pixelSensitivityMask = psf.pixelSensitivityMask;
    par.nDiscretizationBFP = psf.nDiscretizationBFP;

    % Fluorophore
    par.shotNoise =  0; 
    par.reducedExcitation = 0;

    % Microscope setup
    par.wavelength = psf.wavelength;
    par.astigmatism = psf.astigmatism;
    par.objectiveNA = psf.objectiveNA;
    par.objectiveFocalLength = psf.objectiveFocalLength;
    par.refractiveIndices = psf.refractiveIndices;
    par.heightIntermediateLayer = psf.heightIntermediateLayer;

    % Back focal plane
    par.phaseMask = psf.phaseMask;
    par.nDiscretizationBFP = psf.nDiscretizationBFP;

    % Camera
    par.pixelSize = psf.pixelSize;

    % dave apr 2025
    par.nPixels = psf.nPixels;

end