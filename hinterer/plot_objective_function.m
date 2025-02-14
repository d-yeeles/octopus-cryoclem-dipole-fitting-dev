
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
    psfEstimate.position = Length([params(1:2), 0], 'nm');
    psfEstimate.defocus = Length(params(3), 'nm');
    psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
    noiseEstimate = psfEstimate.backgroundNoise;
    nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));

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
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_objective_testing/new_sims_highN/';

% Locate all the images and settings files in the dir
files = dir(frames_dir);
valid_extensions_image = {'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
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
for frame_index = 1:1%1:length(frame_paths)%randperm(length(frame_paths), 50)

    tic; % timing each frame

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

    % Generate throwaway PSF object, later replace the image with our masked image
    par.position = Length([0 0 0], 'nm');
    psfInit = PSF(par);
    psfInit.image = psf_image;
    


    % Example x and y values (fixing them to a specific point)
    x = positionX_nm_array;  % Example value
    y = positionY_nm_array;  % Example value
    defocus = 0;

    % Azimuth values to vary
    azVals = linspace(0, 2*pi, 50);  % Azimuth values, from 0 to 2*pi
    
    % Inclination values (in radians)
    inclinations = [0, 22.5*pi/180, 45*pi/180, 67.5*pi/180, 90*pi/180];
    
    % Create a new figure
    figure;
    
    % Loop over each inclination and calculate the corresponding cost function
    hold on;  % Keep the same plot for multiple curves
    
    for i = 1:length(inclinations)
        inclination = inclinations(i);
        
        % Initialize cost function values for this inclination
        costValues = zeros(length(azVals), 1);
        
        % Calculate the cost function for each azimuth for the current inclination
        for j = 1:length(azVals)
            az = azVals(j);
        
            % Calculate the PSF and cost function
            currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, az]);
        
            % Calculate the cost value (negative log-likelihood)
            costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
            costValues(j) = costValue;  % Store the cost value
        end
        
        % Plot the cost function for this inclination
        plot(azVals*180/pi, costValues, 'LineWidth', 2, 'DisplayName', sprintf('θ = %.2f°', inclination * 180 / pi));
    end
    
    % Customize plot
    xlabel('Azimuth, deg');
    ylabel('Objective function');
    title('Objective function vs Azimuth for Different Inclinations');
    legend show;  % Display legend
    grid on;
    hold off;  % End holding the plot for multiple curves



    
    % % Plot XY + obj func
    % % Plot grid
    % xVals = linspace(positionX_nm_array-50, positionX_nm_array+50, 50);
    % yVals = linspace(positionY_nm_array-50, positionY_nm_array+50, 50);
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % inclination = angleInclination_array;
    % azimuth = angleAzimuth_array;
    % 
    % % Initialize cost function values
    % costMatrix = zeros(length(xVals), length(yVals));
    % 
    % % Compute cost function for each inclination value and place on subplot
    % 
    % tic;
    % 
    % % Initialize cost function values
    % costMatrix = zeros(length(xVals), length(yVals));
    % 
    % % Compute the cost function over x and y coordinates
    % for i = 1:length(xVals)
    %     for j = 1:length(yVals)
    %         % Current (x, y) position
    %         x = xVals(i);
    %         y = yVals(j);
    % 
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, azimuth]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % % Create a 3D surface plot
    % figure;
    % surf(xVals, yVals, costMatrix, 'EdgeColor', 'none');
    % view(3);
    % colorbar;
    % colormap jet;
    % xlabel('x, nm)');
    % ylabel('y, nm)');
    % zlabel('obj func');
    % title(sprintf('obj func, θ = %.2f°)', inclination*180/pi, azimuth*180/pi));
    % hold off;
    % 
    % 
    % 
    % % Plot Xϕ + obj func
    % % Plot grid
    % xVals = linspace(positionX_nm_array-50, positionX_nm_array+50, 50);
    % azVals = linspace(0, 2*pi, 50);  % Round to ensure it's a whole number of points
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % inclination = angleInclination_array;
    % y = positionY_nm_array;
    % 
    % % Initialize cost function values
    % costMatrix = zeros(length(xVals), length(azVals));
    % 
    % % Compute cost function for each inclination value and place on subplot
    % 
    % tic;
    % 
    % % Initialize cost function values
    % costMatrix = zeros(length(xVals), length(azVals));
    % 
    % % Compute the cost function over x and y coordinates
    % for i = 1:length(xVals)
    %     for j = 1:length(azVals)
    %         % Current (x, y) position
    %         x = xVals(i);
    %         az = azVals(j);
    % 
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, az]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % % Create a 3D surface plot
    % figure;
    % surf(xVals, azVals*180/pi, costMatrix, 'EdgeColor', 'none');
    % view(3);
    % colorbar;
    % colormap jet;
    % xlabel('x, nm)');
    % ylabel('az, deg)');
    % zlabel('obj func');
    % title(sprintf('obj func, θ = %.2f°)', inclination*180/pi, azimuth*180/pi));
    % hold off;
    % 
    % 
    % 
    % 
    % % Plot Yϕ + obj func
    % % Plot grid
    % yVals = linspace(positionY_nm_array-50, positionY_nm_array+50, 50);
    % azVals = linspace(0, 2*pi, 50);  % Round to ensure it's a whole number of points
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % inclination = angleInclination_array;
    % x = positionX_nm_array;
    % 
    % % Initialize cost function values
    % costMatrix = zeros(length(yVals), length(azVals));
    % 
    % % Compute cost function for each inclination value and place on subplot
    % 
    % tic;
    % 
    % % Initialize cost function values
    % costMatrix = zeros(length(yVals), length(azVals));
    % 
    % % Compute the cost function over x and y coordinates
    % for i = 1:length(yVals)
    %     for j = 1:length(azVals)
    %         % Current (x, y) position
    %         y = yVals(i);
    %         az = azVals(j);
    % 
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, az]); 
    % 
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Negative log-likelihood
    %     end
    % end
    % 
    % % Create a 3D surface plot
    % figure;
    % surf(yVals, azVals*180/pi, costMatrix, 'EdgeColor', 'none');
    % view(3);
    % colorbar;
    % colormap jet;
    % xlabel('y, nm)');
    % ylabel('az, deg)');
    % zlabel('obj func');
    % title(sprintf('obj func, θ = %.2f°)', inclination*180/pi, azimuth*180/pi));
    % hold off;






    % % Plot Rϕ + obj func
    % % Define radial and azimuthal values for the grid
    % rVals = linspace(0, 100, 10);  % Radial distance
    % azVals = linspace(0, 2*pi, 10);  % Azimuthal angle
    % 
    % % Create meshgrid for r and az
    % [R, AZ] = meshgrid(rVals, azVals);
    % 
    % % Convert from polar to Cartesian coordinates for x, y
    % X = R .* cos(AZ);
    % Y = R .* sin(AZ);
    % 
    % % Fix defocus, inclination, and azimuth for cost function
    % defocus = 0;
    % inclination = angleInclination_array;
    % 
    % % Initialize cost function matrix
    % costMatrix = zeros(length(rVals), length(azVals));
    % 
    % % Compute the cost function over the polar grid
    % tic;
    % 
    % for i = 1:length(rVals)
    %     for j = 1:length(azVals)
    %         % Current (x, y) position
    %         x = X(i, j);
    %         y = Y(i, j);
    %         az = azVals(j);
    % 
    %         % Calculate the PSF and cost function
    %         currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, az]);
    % 
    %         % Calculate the cost value (negative log-likelihood)
    %         costValue = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    %         costMatrix(i, j) = costValue;  % Store the cost value
    %     end
    % end
    % 
    % % Create a 3D surface plot
    % figure;
    % surf(R, AZ*180/pi, costMatrix, 'EdgeColor', 'none');
    % view(3);
    % colorbar;
    % colormap jet;
    % xlabel('r, nm');
    % ylabel('Azimuth, deg');
    % zlabel('Objective function');
    % title(sprintf('Objective function, Inclination = %.2f°', inclination * 180 / pi));
    % 
    % hold off;


    % % Plot XYϕ + obj func contour
    % % Plot grid
    % xVals = linspace(positionX_nm_array-50, positionX_nm_array+50, 10);
    % yVals = linspace(positionY_nm_array-50, positionY_nm_array+50, 10);
    % azVals = linspace(0, 2*pi, 10);  % Round to ensure it's a whole number of points
    % 
    % % Fix defocus, inc, az
    % defocus = 0;
    % inclination = angleInclination_array;
    % azimuth = angleAzimuth_array;
    % 
    % [xGrid, yGrid, azGrid] = meshgrid(xVals, yVals, azVals);
    % costVolume = zeros(size(xGrid));
    % 
    % for i = 1:numel(xGrid)
    %     x = xGrid(i);
    %     y = yGrid(i);
    %     azimuth = azGrid(i);
    %     currentPSF = createFitPSF(psfInit, [x, y, defocus, inclination, azimuth]); 
    %     costVolume(i) = -sum(psf_image .* log(currentPSF) - currentPSF - log(gamma(psf_image + 1)), 'all');
    % end
    % 
    % figure;
    % slice(xGrid, yGrid, azGrid*180/pi, costVolume, [], [], linspace(0, 360, 6));
    % shading interp;
    % colorbar;
    % xlabel('x (nm)');
    % ylabel('y (nm)');
    % zlabel('Azimuth (°)');
    % title('Cost Function Slices in 3D');
    % colormap jet;





    % % Print ground truth obj func value
    % ground_truth_params = [positionX_nm_array, positionY_nm_array, 0, angleInclination_array, angleAzimuth_array];
    % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
    % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');
    % disp(objective_function_true)

    % % Clip values just for display
    % display_image = (currentPSF - min(currentPSF(:))) / (max(currentPSF(:)) - min(currentPSF(:)));
    % imshow(display_image)

    elapsed_time = toc;
    fprintf('    Time to plot: %.2f seconds\n', elapsed_time);


end % end loop over frames
