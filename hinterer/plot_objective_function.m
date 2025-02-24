
%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Fit
%% ----------


% Function to create the current PSF (same as `createFitPSF` method in the class)
function currentPSF = createFitPSF(psfEstimate, params)
    psfEstimate.position = Length([params(1), params(2), 0], 'nm');
    psfEstimate.defocus = Length(params(3), 'nm');
    psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
    noiseEstimate = 3;
    nPhotonEstimate = 1e6;%round(sum(sum(psfEstimate.image - noiseEstimate)));

    % Simulate the PSF
    bfp = BackFocalPlane(psfEstimate);
    % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
    psfEstimate.backFocalPlane = bfp;

    % Apply phase mask
    psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
    psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);

    % Apply attenuation mask
    psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
    psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);

    currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 

    for k=1:size(psfEstimate.stageDrift.motion,1)
        % Apply aberrations
        aberrations = getAberrations(psfEstimate,k);
        aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
        
        % Get image from BFP field
        currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
    end

    currentPsf = adjustExcitation(psfEstimate, currentPsf);
    currentPsf = applyShotNoise(psfEstimate, currentPsf);
    currentPsf = addBackgroundNoise(psfEstimate, currentPsf);

    totalIntensity = sum(currentPsf,'all');
    currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
    currentFitPSF = currentPsf ./ norm(currentPsf);

    currentPSF = currentFitPSF;
end


% This for if checking reparameterised version
% Function to create the current PSF (same as `createFitPSF` method in the class)
function currentPSF = createFitPSF_reparam(psfEstimate, params)


    psfEstimate.position = Length([params(1:2), 0], 'nm');
    psfEstimate.defocus = Length(params(3), 'nm');

    inclination = 0.5*asin(params(4));
    azimuth = atan2(params(6), params(5));

    % disp(lateralPositionAndDefocus(4))

    inclination = mod(inclination, pi/2);
    azimuth = mod(azimuth, 2*pi);

    psfEstimate.dipole = Dipole(inclination, azimuth); % dave jan 2025 - adding angle optimiser

    noiseEstimate = 3;
    nPhotonEstimate = 1e6;%round(sum(sum(psfEstimate.image - noiseEstimate)));

    % Simulate the PSF
    bfp = BackFocalPlane(psfEstimate);
    % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
    psfEstimate.backFocalPlane = bfp;

    % Apply phase mask
    psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
    psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);

    % Apply attenuation mask
    psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
    psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);

    currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 

    for k=1:size(psfEstimate.stageDrift.motion,1)
        % Apply aberrations
        aberrations = getAberrations(psfEstimate,k);
        aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
        
        % Get image from BFP field
        currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
    end

    currentPsf = adjustExcitation(psfEstimate, currentPsf);
    currentPsf = applyShotNoise(psfEstimate, currentPsf);
    currentPsf = addBackgroundNoise(psfEstimate, currentPsf);

    totalIntensity = sum(currentPsf,'all');
    currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
    currentFitPSF = currentPsf ./ norm(currentPsf);

    currentPSF = currentFitPSF;
end




% Input params
% frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_objective_testing/sims/blank/';
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/background_images/highN/';

% Locate all the images and settings files in the dir
files = dir(frames_dir);
valid_extensions_image = {'.tif'};
valid_extensions_settings = {'.m'};
frame_paths = {};
settings_paths = {};

for i = 1:length(files)
    [~, ~, ext] = fileparts(files(i).name);
    if ismember(lower(ext), valid_extensions_image)
        frame_paths{end+1} = fullfile(frames_dir, files(i).name);
    elseif ismember(lower(ext), valid_extensions_settings)
        settings_paths{end+1} = fullfile(frames_dir, files(i).name);
    end
end

% Create the centroids in image coordinates
% 'Read in' ground truth positions (this will be rpelaced by like reading in thunderstorm estimate or something)

% Loop over each image path and process

% Loop over each frame
for frame_index = 4:4%1:length(frame_paths)%randperm(length(frame_paths), 50)

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, length(frame_paths));
    frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};
    
    % Read in ground truth data
    % run(settings_path);
    evalc('run(settings_path)');
    
    par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = image_size_px;
    par.wavelength = Length(wavelength,'nm');
    
    % Load frame
    psf_image = imread(frame_path);
    
    % % Convert to greyscale if not
    % if size(psf_image, 3) == 3
    %     psf_image = rgb2gray(psf_image);
    % end
    
    psf_image = double(psf_image); % Convert to double for calculations
    % psf_image = max(min(psf_image, 255), 0); % Clip values to [0, 255]
    % psf_image = uint8(psf_image); % Convert back to uint8
    % psf_image = double(psf_image);
    % psf_image = 10 + ((psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)))); % need to do this to normalise it to 0,1
    psf_image = psf_image ./ norm(psf_image);
    
    % % Clip values just for display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)
    
    
    
    % % Plot XY, obj
    % 
    % for inclination = [0*pi/180, 22.5*pi/180, 45*pi/180, 67.5*pi/180, 90*pi/180]
    % 
    %     % Generate throwaway PSF object, later replace the image with our masked image
    %     par.position = Length([0, 0, 0], 'nm');
    %     par.dipole = Dipole(inclination, 0);
    %     psfInit = PSF(par);
    %     psfInit.image = psf_image;
    % 
    %     % Plot XY + obj func
    %     % Plot grid
    %     extent = 50;
    %     xVals = linspace(-50, 50, 30);
    %     yVals = linspace(-50, 50, 30);
    % 
    %     % Fix defocus, inc, az
    %     defocus = 0;
    %     azimuth = 0;
    % 
    %     costMatrix = zeros(length(xVals), length(yVals));
    % 
    %     tic;
    % 
    %     for i = 1:length(xVals)
    %         for j = 1:length(yVals)
    % 
    %             x = xVals(i);
    %             y = yVals(j);
    % 
    %             currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, azimuth]); 
    % 
    %             costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %             costMatrix(i, j) = costValue;  % Negative log-likelihood
    %         end
    %     end
    % 
    %     figure;
    %     surf(xVals, yVals, costMatrix, 'EdgeColor', 'none');
    %     view(3);%0, 90);
    %     colorbar;
    %     colormap jet;
    %     xlim([min(xVals) max(xVals)]);
    %     ylim([min(yVals) max(yVals)]);
    %     % zlim([0 300]);
    %     xlabel('x, nm)');
    %     ylabel('y, nm)');
    %     zlabel('obj func');
    %     title(sprintf('obj func, θ = %.2f°)', inclination*180/pi));
    %     hold off;
    % 
    %     % Find minimum value in the costMatrix and its corresponding x and y
    %     [minCostValue, linearIndex] = min(costMatrix(:));
    %     [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    %     % Get the corresponding x and y values
    %     minX = xVals(row);
    %     minY = yVals(col);
    % 
    %     % Display the results
    %     disp(['(', num2str(minX), ' , ', num2str(minY), ')']);
    % 
    %     % % Print ground truth obj func value
    %     % ground_truth_params = [positionX_nm_array, positionY_nm_array, 0, angleInclination_array, angleAzimuth_array];
    %     % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
    %     % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');
    %     % disp(objective_function_true)
    % 
    %     % % Clip values just for display
    %     % display_image = (currentPSF - min(currentPSF(:))) / (max(currentPSF(:)) - min(currentPSF(:)));
    %     % imshow(display_image)
    % 
    %     elapsed_time = toc;
    %     fprintf('    Time to plot: %.2f seconds\n', elapsed_time);
    % 
    % 
    % end % end loop over frames
    
    
    
    % Generate throwaway PSF object, later replace the image with our masked image
    par.position = Length([0, 0, 0], 'nm');
    par.dipole = Dipole(0, 0);
    psfInit = PSF(par);
    psfInit.image = psf_image;
    
    % % Plot Xtheta, obj
    % 
    % % Plot grid
    % extent = 50;
    % xVals = linspace(-20, 20, 320);
    % % incVals = linspace(0, pi/2, 90);
    % incVals = linspace(50*pi/180, 70*pi/180, 320);
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % y = 0;
    % azimuth = 0;
    % 
    % costMatrix = zeros(length(xVals), length(incVals));
    % 
    % for i = 1:length(xVals)
    %     for j = 1:length(incVals)
    % 
    %         x = xVals(i);
    %         inclination = incVals(j);
    % 
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, azimuth]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % figure;
    % surf(xVals, incVals*180/pi, costMatrix', 'EdgeColor', 'none');
    % view(3);%0, 90);
    % colorbar;
    % colormap jet;
    % xlim([min(xVals) max(xVals)]);
    % ylim([min(incVals*180/pi) max(incVals*180/pi)]);
    % % zlim([0 300]);
    % xlabel('x, nm)');
    % ylabel('inc, deg)');
    % zlabel('obj func');
    % title('obj func over x and θ)');
    % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = xVals(row);
    % minY = incVals(col);
    % 
    % % Display the results
    % disp(['(', num2str(minX), ' , ', num2str(minY*180/pi), ')']);
    
    
    
    % % Plot Ytheta, obj
    % 
    % % Plot grid
    % extent = 50;
    % yVals = linspace(-20, 20, 40);
    % incVals = linspace(50*pi/180, 70*pi/180, 20);
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % x = 0;
    % azimuth = 0;
    % 
    % costMatrix = zeros(length(yVals), length(incVals));
    % 
    % for i = 1:length(yVals)
    %     for j = 1:length(incVals)
    % 
    %         y = yVals(i);
    %         inclination = incVals(j);
    % 
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, azimuth]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % figure;
    % surf(yVals, incVals*180/pi, costMatrix', 'EdgeColor', 'none');
    % view(3);%0, 90);
    % colorbar;
    % colormap jet;
    % xlim([min(yVals) max(yVals)]);
    % ylim([min(incVals*180/pi) max(incVals*180/pi)]);
    % % zlim([0 300]);
    % xlabel('y, nm)');
    % ylabel('inc, deg)');
    % zlabel('obj func');
    % title('obj func over y and θ)');
    % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = yVals(row);
    % minY = incVals(col);
    % 
    % % Display the results
    % disp(['(', num2str(minX), ' , ', num2str(minY*180/pi), ')']);
    
    
    
    
    % % Plot theta/phi, obj
    % 
    % % Plot grid
    % incVals = linspace(0, pi/2, 80);
    % azVals = linspace(0, 2*pi, 80);
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % x = 0;
    % y = 0;
    % 
    % costMatrix = zeros(length(incVals), length(azVals));
    % 
    % for i = 1:length(incVals)
    %     for j = 1:length(azVals)
    % 
    %         inclination = incVals(i);
    %         azimuth = azVals(j);
    % 
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, azimuth]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % % figure;
    % % surf(incVals*180/pi, azVals*180/pi, costMatrix', 'EdgeColor', 'none');
    % % view(3);%0, 90);
    % % colorbar;
    % % colormap jet;
    % % xlim([min(incVals*180/pi) max(incVals*180/pi)]);
    % % ylim([min(azVals*180/pi) max(azVals*180/pi)]);
    % % % zlim([0 300]);
    % % xlabel('inc, deg)');
    % % ylabel('az, deg)');
    % % zlabel('obj func');
    % % title('obj func over y and θ)');
    % % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = incVals(row);
    % minY = azVals(col);
    % 
    % % Display the results
    % disp(['(', num2str(minX), ' , ', num2str(minY*180/pi), ')']);
    % 
    % 
    % 
    % 
    % 
    % % Number of patches in the azimuth direction
    % numPatches = 3;  % You can change this to stack more patches
    % 
    % figure;
    % hold on;
    % 
    % for patchIndex = 0:numPatches-1
    %     % Shift the azimuth values for each patch
    %     shiftedAzVals = azVals + patchIndex * 2*pi;  
    % 
    %     % Mirror the inclination values
    %     mirroredIncVals = [incVals, 2*(pi/2) - fliplr(incVals)];
    % 
    %     % Mirror the cost matrix vertically
    %     mirroredCostMatrix = [costMatrix; flipud(costMatrix)];
    % 
    %     % Plot each shifted and mirrored patch
    %     surf(mirroredIncVals * 180/pi, shiftedAzVals * 180/pi, mirroredCostMatrix', 'EdgeColor', 'none');
    % end
    % 
    % view(3);
    % colorbar;
    % colormap jet;
    % xlim([min(incVals*180/pi), 2*max(incVals*180/pi)]);
    % ylim([min(azVals*180/pi), max(azVals*180/pi) + (numPatches-1) * 360]); % Expand range for stacked patches
    % xlabel('Inclination (deg)');
    % ylabel('Azimuth (deg)');
    % zlabel('Objective Function');
    % title('Stacked & Mirrored Cost Function');
    % hold off;





    % But I've now reparameterised to optimise over
    % sin(phi), cos(phi), and sin(2theta)
    % so we need to see how they look really


    % Plot sin(phi)/sin(2*theta), obj

    sinazVals = linspace(-1, 1, 20);
    sinincVals = linspace(-1, 1, 20);

    % Fix defocus, inc, az
    x = 0;
    y = 0;
    defocus = 0;
    cosAz = 1;

    costMatrix = zeros(length(sinazVals), length(sinincVals));

    for i = 1:length(sinazVals)
        for j = 1:length(sinincVals)

            sinAz = sinazVals(i);
            sinInc = sinincVals(j);

            currentPSF = createFitPSF_reparam(psfInit, [x, y, defocus, sinAz, cosAz, sinInc]); 

            costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
            costMatrix(i, j) = costValue;  % Negative log-likelihood
        end
    end

    figure;
    surf(sinazVals, sinincVals, costMatrix', 'EdgeColor', 'none');
    view(3);%0, 90);
    colorbar;
    colormap jet;
    xlim([min(sinazVals) max(sinincVals)]);
    ylim([min(sinincVals) max(sinincVals)]);
    % zlim([0 300]);
    xlabel('sin(az)');
    ylabel('sin(inc)');
    zlabel('obj func');
    title('obj func over sin(az) and sin(inc))');
    hold off;

    % Find minimum value in the costMatrix and its corresponding x and y
    [minCostValue, linearIndex] = min(costMatrix(:));
    [row, col] = ind2sub(size(costMatrix), linearIndex);

    % Get the corresponding x and y values
    minX = sinazVals(row);
    minY = sinincVals(col);

    % Display the results
    disp(['(', num2str((asin(minX))*180/pi), ' , ', num2str((0.5*asin(minY))*180/pi), ')']);









    % % Plot cos(phi)/sin(2*theta), obj
    % 
    % cosazVals = linspace(-1, 1, 20);
    % sinincVals = linspace(-1, 1, 20);
    % 
    % % Fix defocus, inc, az
    % x = 0;
    % y = 0;
    % defocus = 0;
    % sinAz = 1;
    % 
    % costMatrix = zeros(length(cosazVals), length(sinincVals));
    % 
    % for i = 1:length(cosazVals)
    %     for j = 1:length(sinincVals)
    % 
    %         cosAz = cosazVals(i);
    %         sinInc = sinincVals(j);
    % 
    %         currentPSF = createFitPSF_reparam(psfInit, [x, y, defocus, sinAz, cosAz, sinInc]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % figure;
    % surf(cosazVals, sinincVals, costMatrix', 'EdgeColor', 'none');
    % view(3);%0, 90);
    % colorbar;
    % colormap jet;
    % xlim([min(cosazVals) max(sinincVals)]);
    % ylim([min(cosazVals) max(sinincVals)]);
    % % zlim([0 300]);
    % xlabel('cos(az)');
    % ylabel('sin(inc)');
    % zlabel('obj func');
    % title('obj func over cos(az) and sin(inc))');
    % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = cosazVals(row);
    % minY = sinincVals(col);
    % 
    % % Display the results
    % disp(['(', num2str((acos(minX))*180/pi), ' , ', num2str((0.5*asin(minY))*180/pi), ')']);
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % % Plot cos(phi)/sin(2*theta), obj
    % 
    % cosazVals = linspace(-1, 1, 20);
    % sinazVals = linspace(-1, 1, 20);
    % 
    % % Fix defocus, inc, az
    % x = 0;
    % y = 0;
    % defocus = 0;
    % sinInc = 1;
    % 
    % costMatrix = zeros(length(cosazVals), length(sinazVals));
    % 
    % for i = 1:length(cosazVals)
    %     for j = 1:length(sinazVals)
    % 
    %         cosAz = cosazVals(i);
    %         sinAz = sinazVals(j);
    % 
    %         currentPSF = createFitPSF_reparam(psfInit, [x, y, defocus, sinAz, cosAz, sinInc]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % figure;
    % surf(cosazVals, sinazVals, costMatrix', 'EdgeColor', 'none');
    % view(3);%0, 90);
    % colorbar;
    % colormap jet;
    % xlim([min(cosazVals) max(sinazVals)]);
    % ylim([min(cosazVals) max(sinincVals)]);
    % % zlim([0 300]);
    % xlabel('cos(az)');
    % ylabel('sin(az)');
    % zlabel('obj func');
    % title('obj func over cos(az) and sin(inc))');
    % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = cosazVals(row);
    % minY = sinazVals(col);
    % 
    % % Display the results
    % disp(['(', num2str((acos(minX))*180/pi), ' , ', num2str((0.5*asin(minY))*180/pi), ')']);




    % 
    % 
    % % Plot az/sin(2*theta), obj
    % 
    % azVals = linspace(0, 2*pi, 20);
    % sinincVals = linspace(-1, 1, 20);
    % 
    % % Fix defocus, inc, az
    % x = 0;
    % y = 0;
    % defocus = 0;
    % sinInc = 1;
    % 
    % costMatrix = zeros(length(azVals), length(sinincVals));
    % 
    % for i = 1:length(azVals)
    %     for j = 1:length(sinincVals)
    % 
    %         azimuth = azVals(i);
    %         sinInc = sinincVals(j);
    % 
    %         currentPSF = createFitPSF_reparam(psfInit, [x, y, defocus, sinInc, azimuth]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % figure;
    % surf(azVals*180/pi, sinincVals, costMatrix', 'EdgeColor', 'none');
    % view(3);%0, 90);
    % colorbar;
    % colormap jet;
    % xlim([min(azVals*180/pi) max(azVals*180/pi)]);
    % ylim([min(sinincVals) max(sinincVals)]);
    % % zlim([0 300]);
    % xlabel('az');
    % ylabel('sin(az)');
    % zlabel('obj func');
    % title('obj func over az and sin(inc))');
    % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = azVals(row);
    % minY = sinincVals(col);
    % 
    % % Display the results
    % disp(['(', num2str(minX*180/pi), ' , ', num2str((0.5*asin(minY))*180/pi), ')']);
    % 




    % 
    % % Plot az/sin(2*theta), obj
    % 
    % incVals = linspace(0, 2*pi, 20);
    % sinazVals = linspace(-1, 1, 20);
    % 
    % % Fix defocus, inc, az
    % x = 0;
    % y = 0;
    % defocus = 0;
    % cosAz = 1;
    % 
    % costMatrix = zeros(length(incVals), length(sinazVals));
    % 
    % for i = 1:length(incVals)
    %     for j = 1:length(sinazVals)
    % 
    %         inclination = incVals(i);
    %         sinAz = sinazVals(j);
    % 
    %         currentPSF = createFitPSF_reparam(psfInit, [x, y, defocus, inclination, cosAz, sinAz]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % figure;
    % surf(incVals*180/pi, sinazVals, costMatrix', 'EdgeColor', 'none');
    % view(3);%0, 90);
    % colorbar;
    % colormap jet;
    % xlim([min(incVals*180/pi) max(incVals*180/pi)]);
    % ylim([min(sinazVals) max(sinazVals)]);
    % % zlim([0 300]);
    % xlabel('inc');
    % ylabel('sin(az)');
    % zlabel('obj func');
    % title('obj func over az and sin(az))');
    % hold off;
    % 
    % % Find minimum value in the costMatrix and its corresponding x and y
    % [minCostValue, linearIndex] = min(costMatrix(:));
    % [row, col] = ind2sub(size(costMatrix), linearIndex);
    % 
    % % Get the corresponding x and y values
    % minX = incVals(row);
    % minY = sinazVals(col);
    % 
    % % Display the results
    % disp(['(', num2str(minX*180/pi), ' , ', num2str((asin(minY))*180/pi), ')']);






end