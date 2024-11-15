%% Example script for simulating and fitting PSF

close all
clear all

% Create a tiled layout
t = tiledlayout(4, 4);  % 2 rows, 2 columns

% Loop to create plots
for i = 0:15  % Loop from 0 to 8
    nexttile;  % Move to the next tile
    
    % Call the function to generate and plot the PSF
    plotPSF(i);
end

% Function to plot PSF
function plotPSF(i)
    % Set input parameters
    angleInclination = i * pi / 16;
    angleAzimuth = i * pi / 8;
    par.position = Length([i 0 0] * 100, 'nm');
    par.dipole = Dipole(angleInclination, angleAzimuth);
    par.stageDrift = LinearDrift(0.0001, 0.0001, 1, 'nm')%(10, 1, 10, 'nm');
    par.defocus = Length(-500 + 1000 * rand(), 'nm');
    
    % Simulate psf
    psf = PSF(par);
    
    % Fitting
    parEst.angleInclinationEstimate = angleInclination;
    parEst.angleAzimuthEstimate = angleAzimuth;
    parEst.stageDrift = par.stageDrift;
    fitResult = FitPSF(psf, parEst);
    
    % Plotting the results
    plot(fitResult);
    title(['AngleI: ', num2str(angleInclination), 'AngleA: ', num2str(angleAzimuth)]);
end

% Add a common title for the tiled layout
title(t, 'PSF fitting on dipoles with varying inc, az');



