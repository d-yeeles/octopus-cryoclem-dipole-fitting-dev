%% Example script for simulating and fitting PSF

close all
clear all

addpath(genpath('../'));


% Set input parameters
angleInclination = 0;
x_true = 900;
y_true = 900;
pixel_size_nm = 51.2;
image_size_px = 51; % must be odd



par.position = Length([x_true y_true 0],'nm');
par.dipole = Dipole(angleInclination, 0);
% par.stageDrift = LinearDrift(10,1,10,'nm');
par.defocus = Length(-500+1000*rand(), 'nm');
par.pixelSize = Length(pixel_size_nm,'nm');
par.nPixels = image_size_px;

% attocube params
par.wavelength = Length(500,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(770,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
% par.backgroundNoise = 10; % taken from looking at blank bit of example data
par.nPhotons = 1e10;%1e20;

% Simulate psf
tic; % Start the timer
psf = PSF(par);
toc; % Stop the timer and display the elapsed time

% Fitting
parEst.angleInclinationEstimate = angleInclination;
parEst.angleAzimuthEstimate = 0;
% parEst.stageDrift = par.stageDrift;
tic; % Start the timer
fitResult = FitPSF(psf, parEst);
toc; % Stop the timer and display the elapsed time

% Print localisations
x_guess = fitResult.estimatesPositionDefocus.ML(1);
y_guess = fitResult.estimatesPositionDefocus.ML(2);

disp(['true  = (', num2str(x_true), ', ', num2str(y_true), ')'])
disp(['guess = (', num2str(x_guess), ', ', num2str(y_guess), ')'])

% Plot
figure
plot(fitResult)
addMM=@(x) sprintf('%i',x*pixel_size_nm - image_size_px*pixel_size_nm/2);
xticklabels(cellfun(addMM,num2cell(xticks'),'UniformOutput',false));
yticklabels(cellfun(addMM,num2cell(yticks'),'UniformOutput',false));
% figure
% imagesc(fitResult.image)
% hold on;
% plot(x_guess, y_guess, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
% plot(x_true, y_true, 'x', 'MarkerFaceColor', 'r', 'MarkerSize', 10);