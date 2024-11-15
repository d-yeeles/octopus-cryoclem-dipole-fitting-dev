%% Example script for simulating and fitting multiple PSFs

close all
clear all

addpath(genpath('./'));

% Define base parameters
baseParams.nPixels = 171;
baseParams.stageDrift = LinearDrift(0.0001, 0.0001, 1, 'nm');
baseParams.defocusRange = [-500, 500];  % Defocus range

% Define angles and positions for each PSF you want to create
angleIncs = [pi/16, pi/8];  % Inclination angles for each PSF
positions = [0, 50; 0, 0; 0, 0];  % Positions for each PSF (each column is a PSF)

% Create a figure to display the PSFs
figure;
hold on;  % Hold on to overlay multiple plots

% Loop to generate and plot each PSF
for i = 1:length(angleIncs)
    % Set specific parameters for this PSF
    par = baseParams;  % Start with base parameters
    par.position = Length(positions(:, i)' .* 100, 'nm');
    par.dipole = Dipole(angleIncs(i), 0);
    par.defocus = Length(baseParams.defocusRange(1) + diff(baseParams.defocusRange) * rand(), 'nm');

    % Simulate PSF
    psf(i) = PSF(par);

    % Plot with transparency to combine
    imagesc(psf(i).image, 'AlphaData', 0.5);  % Adjust AlphaData for transparency
end

% Add plot details
colorbar;  % Add color bar for reference
title('Combined Heat Maps');  % Title for the combined image
axis equal tight;  % Equal scaling and tight axes
xlabel('X-axis');  % Replace with your actual label
ylabel('Y-axis');

hold off;  % Release hold
