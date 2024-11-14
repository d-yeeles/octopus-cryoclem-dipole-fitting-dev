%% Fitting multiple PSFs in a single frame

close all
clear all

addpath(genpath('./'));

% Define base parameters
%baseParams.nPixels = 171;
%baseParams.stageDrift = LinearDrift(0.0001, 0.0001, 1, 'nm');
%baseParams.defocusRange = [-500, 500];  % Defocus range

% Loop over each blob in a frame
for i = 1:10%length(blobs)

    % --------------------
    % Simulate psf - replace this whole chunk with reading in an image
    % --------------------
    % Not sure what FitPSF() uses to do the fitting. Does it just use the
    % psf.image part of a psf object? Or does it use the other stuff stored in that?
    % Looking at FitPSF.m, it seems to use a lot of stuff that isn't just
    % the image. Like psf.pixelSensitivityMask, psf.stageDrift, and the
    % parameter estimates. So might need to adapt the PSF simulator so that
    % it "simulates" a psf object containing the image we choose to give
    % it, which will be a patch of a frame.
    %
    % Maybe we can adapt the PSF() class/function/whatever that currently
    % generates a simulation so that it simply generates a psf object where
    % the psf.image bit = some input image. The other stuff in the psf
    % object are experimental params mostly? I think.
    %
    % If we use thunderSTORM to get an initial etimate of localisations ,we
    % can pass those as our initial estimates for this to use.
    % Oh and because the patches we consider will be centered on those
    % locations, the initial estimate will always be (0,0) by definition.

    % Set input parameters
    angleInclination = i * pi / 16;
    angleAzimuth = i * pi / 8;
    par.position = Length([i 0 0] * 100, 'nm');
    par.dipole = Dipole(angleInclination, angleAzimuth);
    par.stageDrift = LinearDrift(0.0001, 0.0001, 1, 'nm');%(10, 1, 10, 'nm');
    par.defocus = Length(-500 + 1000 * rand(), 'nm');
    % Simulate
    psf = PSF(par);
    % --------------------

    % Fitting
    parEst.angleInclinationEstimate = angleInclination;
    parEst.angleAzimuthEstimate = angleAzimuth;
    parEst.stageDrift = par.stageDrift;
    fitResult = FitPSF(psf, parEst);

    % Store results

    % Estimates for angles are stored in:
    %   fitResult.angleInclinationEstimate
    %   fitResult.angleAzimuthEstimate
    %
    % Position is estimated with both least squares and max likelihood
    %   fitResult.estimatesPositionDefocus.LS
    %   fitResult.estimatesPositionDefocus.ML
    %       given as [x, y, defocus]

    inclination_ests(i) = fitResult.angleInclinationEstimate;
    azimuth_ests(i) = fitResult.angleAzimuthEstimate;
    x_ests(i) = fitResult.estimatesPositionDefocus.ML(1);
    y_ests(i) = fitResult.estimatesPositionDefocus.ML(2);
    defocus_ests(i) = fitResult.estimatesPositionDefocus.ML(3);

end

disp(['Inclination angles: ', num2str(inclination_ests)]);
disp(['Azimuthal angles: ', num2str(azimuth_ests)]);
disp(['x position: ', num2str(x_ests)]);
disp(['y position: ', num2str(y_ests)]);
disp(['defocus: ', num2str(defocus_ests)]);