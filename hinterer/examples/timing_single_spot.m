%% Example script for simulating and fitting PSF

close all
clear all

addpath(genpath('./'));

% Set input parameters
angleInclination = pi/16;
par.position = Length([0 0 0].*100,'nm');
par.dipole = Dipole(angleInclination, 0);
par.stageDrift = LinearDrift(10,1,10,'nm');
par.defocus = Length(-500+1000*rand(), 'nm');

% Simulate psf
numIterations = 1000;
elapsedTimes = zeros(1, numIterations);
for i = 1:numIterations
    tic;
    psf = PSF(par);
    elapsedTimes(i) = toc;
end
fprintf('Average execution time: %.4f seconds\n', mean(elapsedTimes));

% Fitting
parEst.angleInclinationEstimate = angleInclination;
parEst.angleAzimuthEstimate = 0;
parEst.stageDrift = par.stageDrift;
numIterations = 100;
elapsedTimes = zeros(1, numIterations);
for i = 1:numIterations
    tic;
    fitResult = FitPSF(psf, parEst);
    elapsedTimes(i) = toc;
end
fprintf('Average execution time: %.4f seconds\n', mean(elapsedTimes));
