%% ----------------------------------------------------------------------
%% HoM
%% ----------------------------------------------------------------------

%% Thunderstorm patch version


%% Fitting multiple PSFs in a single frame

% Same as loop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate thpem.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_all/";
stack_path = stack_dir + "sim_stack.tif";

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/fitting_results_hinterer_on_mortensen_all.py';

% Model to use
model = 'hinterer';

% Hard-coded microscope parameters
pixel_size_nm = 52;
wavelength = 500;
objectiveNA = 2.17;
objectiveFocalLength = 770;
refractiveIndices = [1.31 , 2.17 , 2.17];
backgroundNoise = 0;
nPhotons = 2000;
pixelSensitivityMask = PixelSensitivity.uniform(9);
nDiscretizationBFP = 129;

patch_width_nm = 988; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
image_width_px = stack_info(1).Width;
image_height_px = stack_info(1).Height;
number_of_frames = numel(stack_info);

% Function to create the current PSF (same as `createFitPSF` method in the class)
% function currentPSF = createFitPSF(psfEstimate, params)
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % Simulate the PSF
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % Apply phase mask
%     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
%     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% 
%     % Apply attenuation mask
%     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
%     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         % Apply aberrations
%         aberrations = getAberrations(psfEstimate,k);
%         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% 
%         % Get image from BFP field
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
%     end
% 
%     currentPsf = adjustExcitation(psfEstimate, currentPsf);
%     currentPsf = applyShotNoise(psfEstimate, currentPsf);
%     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% end
% function currentPSF = createFitPSF(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end
% 
% function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % bfp = BackFocalPlane(psfEstimate);
%     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end




% % Locate all the images and settings files in the dir
% files = dir(frames_dir);
% valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
% valid_extensions_settings = {'.m'};
% frame_paths = {};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_image)
%         frame_paths{end+1} = fullfile(frames_dir, files(i).name);
%     elseif ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(frames_dir, files(i).name);
%     end
% end

% Locate all the settings files in the dir
% files = dir(stack_dir);
% valid_extensions_settings = {'.m'};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(stack_dir, files(i).name);
%     end
% end

files = dir(fullfile(stack_dir, '*.m')); % Use pattern matching
settings_paths = {};
for i = 1:length(files)
    settings_paths{end+1} = fullfile(stack_dir, files(i).name);
end

% Sort the paths based on frame numbers
frame_numbers = zeros(length(settings_paths), 1);
for i = 1:length(settings_paths)
    [~, name, ~] = fileparts(settings_paths{i});
    % Extract the number after 'params_frame'
    frame_number = str2double(regexp(name, 'params_frame(\d+)', 'tokens', 'once'));
    % frame_number = str2double(regexp(name, 'params_(\d+)', 'tokens', 'once'));
    frame_numbers(i) = frame_number;
end

% Sort paths based on frame numbers
[~, idx] = sort(frame_numbers);
settings_paths = settings_paths(idx);

% % Read in thunderstorm localisations
% thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');
% 
% % Extract the required columns
% thunderstorm_frame_array = thunderstorm_results.("frame");
% thunderstorm_x_array = thunderstorm_results.("x [nm]");
% thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % Initialise overall results arrays
% positionX_trues = [];
% positionY_trues = [];
% angleInclination_trues = [];
% angleAzimuth_trues = [];
% 
% positionX_estimates = [];
% positionY_estimates = [];
% angleInclination_estimates = [];
% angleAzimuth_estimates = [];
% 
% positionX_errors = [];
% positionY_errors = [];
% angleInclination_errors = [];
% angleAzimuth_errors = [];
% 
% % objective_function_trues = [];
% % objective_function_estimates = [];
% % objective_function_errors = [];

% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in parameters/settings if using simulated to get ground truth
    evalc('run(settings_path)');
    % par.pixelSize = Length(pixel_size_nm,'nm');
    % par.objectiveNA = objectiveNA;
    % par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    % par.nPixels = min(image_width_px, image_height_px);
    % par.wavelength = Length(wavelength,'nm');
    % 
    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    % par.nDiscretizationBFP = 129;
    % patch_width_nm = 1000; % size of patch around blob to consider
    % patch_width_px = patch_width_nm/pixel_size_nm;
    % patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    % par.nPixels = min(image_width_px, image_height_px);
    % par.nPixels = 2*floor(min(image_width_px, image_height_px)/2) + 1; % This enforces odd number of pixels, as required by Hinterer
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Take thunderstorm results for that frame
    % % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    % current_frame_mask = thunderstorm_frame_array == frame_index;
    % current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    % current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(positionX_nm_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(positionX_nm_array));

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        % adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        % adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);
        % 
        % % Calculate patch boundaries using the adjusted center
        % patch_start_x_px = adjusted_centre_x_px - half_width;
        % patch_start_y_px = adjusted_centre_y_px - half_width;
        % patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        % patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        % 
        % % Extract patch from image
        % patch_indices_x = patch_start_x_px:patch_end_x_px;
        % patch_indices_y = patch_start_y_px:patch_end_y_px;
        % patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        % 
        % % Get actual patch dimensions after boundary checking
        % actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        % actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        % 
        % % Calculate ACTUAL center of patch in pixel coordinates
        % actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        % actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        % 
        % % Convert ACTUAL patch center to physical coordinates
        % actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        % actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        % 
        % % Calculate expected position relative to patch center - beause of
        % % the boundary checking stuff, it might not quite be zero anymore
        % 
        % expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        % expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;





        % Better approach - avoid rounding and maintain sub-pixel precision
        % First calculate the integer patch boundaries that ensure the patch stays within image
        patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
        patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
        patch_end_x_px = min(image_width, patch_start_x_px + patch_width_px - 1);
        patch_end_y_px = min(image_height, patch_start_y_px + patch_width_px - 1);
        
        % Adjust patch start if necessary to maintain constant patch size
        if (patch_end_x_px - patch_start_x_px + 1) < patch_width_px
            patch_start_x_px = max(1, patch_end_x_px - patch_width_px + 1);
        end
        if (patch_end_y_px - patch_start_y_px + 1) < patch_width_px
            patch_start_y_px = max(1, patch_end_y_px - patch_width_px + 1);
        end
        
        % Extract patch from image
        patch_indices_x = patch_start_x_px:patch_end_x_px;
        patch_indices_y = patch_start_y_px:patch_end_y_px;
        patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        
        % Get actual patch dimensions after boundary checking
        actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        
        % Calculate ACTUAL center of patch in pixel coordinates
        actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        
        % Convert ACTUAL patch center to physical coordinates
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        
        % Calculate expected position relative to patch center
        expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;
        





        % % DEBUG output to see coordinates
        % fprintf('Initial patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     patch_centre_x_nm, patch_centre_y_nm, patch_centre_x_px, patch_centre_y_px);
        % fprintf('Actual patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     actual_patch_center_x_nm, actual_patch_center_y_nm, actual_patch_center_x_px, actual_patch_center_y_px);

        % % Clip values just if want to display
        % display_image = (patch_psf_image - min(patch_psf_image(:))) / (max(patch_psf_image(:)) - min(patch_psf_image(:)));
        % imshow(display_image)

        % Setup PSF object for fitting
        par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        if strcmpi(model, 'hinterer')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
            psfInit = PSF(par);
        elseif strcmpi(model, 'mortensen')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_mortensen(par);
        elseif strcmpi(model, 'gaussian')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_gaussian(par);
        else
            error('Unknown model type: %s', model);
        end
        
        psfInit.image = patch_psf_image;
        psfInit.nPixels = length(patch_psf_image);

        % % If not masking or cropping, use this
        % % Generate throwaway PSF object, replace the image with our image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % psfInit.image = psf_image;

        % Run the fit - given nothing, fmincon - reparameterised

        % Give it a background noise estimate (let's say you know it to within 10%)
        if psfInit.backgroundNoise > 1e-5
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly. Will need to change this becaue won't know this for real data
        else
            parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
        end

        parEst.parameterBounds.x = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.newangle1 = [-1, 1];%[1/sqrt(2) - 0.0000000001, 1/sqrt(2) + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];%[0 - 0.0000000001, 0 + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give patch centre as initial value. Add some fake noise to
        % simulate imperfect feature detection.
        parEst.parameterStartValues.x = Length(expected_offset_x + (0 + 200 * (rand() - 0.5)), 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y + (0 + 200 * (rand() - 0.5)), 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
        objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of things like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_true = positionX_nm_array(blob_index);
        positionY_nm_true = positionY_nm_array(blob_index);
        angleInclination_true = angleInclination_array(blob_index);
        angleAzimuth_true = angleAzimuth_array(blob_index);
        photons_true = photon_estimate;

        positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
        angleInclination_error = angleInclination_true - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
        photons_error = photons_true - photons_fit_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_trues_frame(blob_index) = positionX_nm_true;
        positionY_nm_trues_frame(blob_index) = positionY_nm_true;
        angleInclination_trues_frame(blob_index) = angleInclination_true;
        angleAzimuth_trues_frame(blob_index) = angleAzimuth_true;
        photons_trues_frame(blob_index) = photons_true;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_fit_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;




        
        
        % % This will create and save a plot for each patch showing the true and estimated positions
        % 
        % % Create figure
        % figure('Position', [100, 100, 600, 500]);
        % 
        % % Display the patch image
        % imagesc(patch_psf_image);
        % colormap('hot'); % Use a colormap that makes the PSF visible
        % axis equal tight;
        % title(sprintf('Frame %d, Blob %d', frame_index, blob_index));
        % colorbar;
        % 
        % % Convert the estimated position from local nm to local pixel coordinates
        % % Note: patch is centered at (half_width+1, half_width+1) in local coordinates
        % local_est_x_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(1) / pixel_size_nm;
        % local_est_y_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(2) / pixel_size_nm;
        % 
        % % Ground truth should be at the center of the patch (approximately)
        % local_gt_x_px = half_width + 1;  % Center of patch
        % local_gt_y_px = half_width + 1;  % Center of patch
        % 
        % hold on;
        % 
        % % Plot true position with a green cross
        % plot(local_gt_x_px, local_gt_y_px, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        % 
        % % Plot estimated position with a red circle
        % plot(local_est_x_px, local_est_y_px, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        % 
        % % Add legend
        % legend('Ground Truth', 'Estimate', 'Location', 'Best');
        % 
        % % Save the figure to the run directory
        % output_filename = sprintf('psf_patch_f%d_b%d.png', frame_index, blob_index);
        % saveas(gcf, fullfile(stack_dir, output_filename));
        % fprintf('Saved visualization to: %s\n', fullfile(stack_dir, output_filename));
        % 
        % % Close the figure
        % close;






        % % DEBUGGING
        % 
        % fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        % fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_true, positionY_nm_true);
        % 
        % % Convert ground truth to pixels
        % gt_x_px = nm_to_px(positionX_nm_true, pixel_size_nm, image_width_px, 'x');
        % gt_y_px = nm_to_px(positionY_nm_true, pixel_size_nm, image_height_px, 'y');
        % fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);
        % 
        % % Show patch info
        % fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
        %         patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        % fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);
        % 
        % % Show actual patch center
        % fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        % fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);
        % 
        % % Calculate ground truth relative to patch center (this should be the target for fitting)
        % fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_true - actual_patch_center_x_nm, ...
        %         positionY_nm_true - actual_patch_center_y_nm);
        % 
        % % After fitting, show the results
        % fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
        %         fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        % fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_estimate, positionY_nm_estimate);
        % fprintf('Error (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_error, positionY_nm_error);
        % fprintf('Error (px): (%.2f, %.2f)\n', ...
        %         positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_trues{frame_index} = positionX_nm_trues_frame;
    positionY_nm_trues{frame_index} = positionY_nm_trues_frame;
    angleInclination_trues{frame_index} = angleInclination_trues_frame;
    angleAzimuth_trues{frame_index} = angleAzimuth_trues_frame;
    photons_trues{frame_index} = photons_trues_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_trues_flat = [positionX_nm_trues{:}];
    positionY_nm_trues_flat = [positionY_nm_trues{:}];
    angleInclination_trues_flat = [angleInclination_trues{:}];
    angleAzimuth_trues_flat = [angleAzimuth_trues{:}];
    photons_trues_flat = [photons_trues{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    photons_errors_flat = [photons_errors{:}];

    % objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_trues_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_trues_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_trues_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_trues_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_trues_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames








%% ----------------------------------------------------------------------
%% MoM
%% ----------------------------------------------------------------------









close all;
clear all;


%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_all/";
stack_path = stack_dir + "sim_stack.tif";

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/fitting_results_mortensen_on_mortensen_all.py';

% Model to use
model = 'mortensen';

% Hard-coded microscope parameters
pixel_size_nm = 52;
wavelength = 500;
objectiveNA = 2.17;
objectiveFocalLength = 770;
refractiveIndices = [1.31 , 2.17 , 2.17];
backgroundNoise = 0;
nPhotons = 2000;
pixelSensitivityMask = PixelSensitivity.uniform(9);
nDiscretizationBFP = 129;

patch_width_nm = 988; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
image_width_px = stack_info(1).Width;
image_height_px = stack_info(1).Height;
number_of_frames = numel(stack_info);

% Function to create the current PSF (same as `createFitPSF` method in the class)
% function currentPSF = createFitPSF(psfEstimate, params)
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % Simulate the PSF
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % Apply phase mask
%     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
%     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% 
%     % Apply attenuation mask
%     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
%     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         % Apply aberrations
%         aberrations = getAberrations(psfEstimate,k);
%         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% 
%         % Get image from BFP field
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
%     end
% 
%     currentPsf = adjustExcitation(psfEstimate, currentPsf);
%     currentPsf = applyShotNoise(psfEstimate, currentPsf);
%     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% end
% function currentPSF = createFitPSF(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end
% 
% function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % bfp = BackFocalPlane(psfEstimate);
%     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end




% % Locate all the images and settings files in the dir
% files = dir(frames_dir);
% valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
% valid_extensions_settings = {'.m'};
% frame_paths = {};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_image)
%         frame_paths{end+1} = fullfile(frames_dir, files(i).name);
%     elseif ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(frames_dir, files(i).name);
%     end
% end

% Locate all the settings files in the dir
% files = dir(stack_dir);
% valid_extensions_settings = {'.m'};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(stack_dir, files(i).name);
%     end
% end

files = dir(fullfile(stack_dir, '*.m')); % Use pattern matching
settings_paths = {};
for i = 1:length(files)
    settings_paths{end+1} = fullfile(stack_dir, files(i).name);
end

% Sort the paths based on frame numbers
frame_numbers = zeros(length(settings_paths), 1);
for i = 1:length(settings_paths)
    [~, name, ~] = fileparts(settings_paths{i});
    % Extract the number after 'params_frame'
    frame_number = str2double(regexp(name, 'params_frame(\d+)', 'tokens', 'once'));
    % frame_number = str2double(regexp(name, 'params_(\d+)', 'tokens', 'once'));
    frame_numbers(i) = frame_number;
end

% Sort paths based on frame numbers
[~, idx] = sort(frame_numbers);
settings_paths = settings_paths(idx);

% % Read in thunderstorm localisations
% thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');
% 
% % Extract the required columns
% thunderstorm_frame_array = thunderstorm_results.("frame");
% thunderstorm_x_array = thunderstorm_results.("x [nm]");
% thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % Initialise overall results arrays
% positionX_trues = [];
% positionY_trues = [];
% angleInclination_trues = [];
% angleAzimuth_trues = [];
% 
% positionX_estimates = [];
% positionY_estimates = [];
% angleInclination_estimates = [];
% angleAzimuth_estimates = [];
% 
% positionX_errors = [];
% positionY_errors = [];
% angleInclination_errors = [];
% angleAzimuth_errors = [];
% 
% % objective_function_trues = [];
% % objective_function_estimates = [];
% % objective_function_errors = [];

% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in parameters/settings if using simulated to get ground truth
    evalc('run(settings_path)');
    % par.pixelSize = Length(pixel_size_nm,'nm');
    % par.objectiveNA = objectiveNA;
    % par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    % par.nPixels = min(image_width_px, image_height_px);
    % par.wavelength = Length(wavelength,'nm');
    % 
    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    % par.nDiscretizationBFP = 129;
    % patch_width_nm = 1000; % size of patch around blob to consider
    % patch_width_px = patch_width_nm/pixel_size_nm;
    % patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    % par.nPixels = min(image_width_px, image_height_px);
    % par.nPixels = 2*floor(min(image_width_px, image_height_px)/2) + 1; % This enforces odd number of pixels, as required by Hinterer
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Take thunderstorm results for that frame
    % % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    % current_frame_mask = thunderstorm_frame_array == frame_index;
    % current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    % current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(positionX_nm_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(positionX_nm_array));

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        % adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        % adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);
        % 
        % % Calculate patch boundaries using the adjusted center
        % patch_start_x_px = adjusted_centre_x_px - half_width;
        % patch_start_y_px = adjusted_centre_y_px - half_width;
        % patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        % patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        % 
        % % Extract patch from image
        % patch_indices_x = patch_start_x_px:patch_end_x_px;
        % patch_indices_y = patch_start_y_px:patch_end_y_px;
        % patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        % 
        % % Get actual patch dimensions after boundary checking
        % actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        % actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        % 
        % % Calculate ACTUAL center of patch in pixel coordinates
        % actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        % actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        % 
        % % Convert ACTUAL patch center to physical coordinates
        % actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        % actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        % 
        % % Calculate expected position relative to patch center - beause of
        % % the boundary checking stuff, it might not quite be zero anymore
        % 
        % expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        % expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;





        % Better approach - avoid rounding and maintain sub-pixel precision
        % First calculate the integer patch boundaries that ensure the patch stays within image
        patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
        patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
        patch_end_x_px = min(image_width, patch_start_x_px + patch_width_px - 1);
        patch_end_y_px = min(image_height, patch_start_y_px + patch_width_px - 1);
        
        % Adjust patch start if necessary to maintain constant patch size
        if (patch_end_x_px - patch_start_x_px + 1) < patch_width_px
            patch_start_x_px = max(1, patch_end_x_px - patch_width_px + 1);
        end
        if (patch_end_y_px - patch_start_y_px + 1) < patch_width_px
            patch_start_y_px = max(1, patch_end_y_px - patch_width_px + 1);
        end
        
        % Extract patch from image
        patch_indices_x = patch_start_x_px:patch_end_x_px;
        patch_indices_y = patch_start_y_px:patch_end_y_px;
        patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        
        % Get actual patch dimensions after boundary checking
        actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        
        % Calculate ACTUAL center of patch in pixel coordinates
        actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        
        % Convert ACTUAL patch center to physical coordinates
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        
        % Calculate expected position relative to patch center
        expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;
        





        % % DEBUG output to see coordinates
        % fprintf('Initial patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     patch_centre_x_nm, patch_centre_y_nm, patch_centre_x_px, patch_centre_y_px);
        % fprintf('Actual patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     actual_patch_center_x_nm, actual_patch_center_y_nm, actual_patch_center_x_px, actual_patch_center_y_px);

        % % Clip values just if want to display
        % display_image = (patch_psf_image - min(patch_psf_image(:))) / (max(patch_psf_image(:)) - min(patch_psf_image(:)));
        % imshow(display_image)

        % Setup PSF object for fitting
        par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        if strcmpi(model, 'hinterer')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
            psfInit = PSF(par);
        elseif strcmpi(model, 'mortensen')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_mortensen(par);
        elseif strcmpi(model, 'gaussian')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_gaussian(par);
        else
            error('Unknown model type: %s', model);
        end
        
        psfInit.image = patch_psf_image;
        psfInit.nPixels = length(patch_psf_image);

        % % If not masking or cropping, use this
        % % Generate throwaway PSF object, replace the image with our image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % psfInit.image = psf_image;

        % Run the fit - given nothing, fmincon - reparameterised

        % Give it a background noise estimate (let's say you know it to within 10%)
        if psfInit.backgroundNoise > 1e-5
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly. Will need to change this becaue won't know this for real data
        else
            parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
        end

        parEst.parameterBounds.x = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.newangle1 = [-1, 1];%[1/sqrt(2) - 0.0000000001, 1/sqrt(2) + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];%[0 - 0.0000000001, 0 + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give patch centre as initial value. Add some fake noise to
        % simulate imperfect feature detection.
        parEst.parameterStartValues.x = Length(expected_offset_x + (0 + 200 * (rand() - 0.5)), 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y + (0 + 200 * (rand() - 0.5)), 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
        objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of things like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_true = positionX_nm_array(blob_index);
        positionY_nm_true = positionY_nm_array(blob_index);
        angleInclination_true = angleInclination_array(blob_index);
        angleAzimuth_true = angleAzimuth_array(blob_index);
        photons_true = photon_estimate;

        positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
        angleInclination_error = angleInclination_true - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
        photons_error = photons_true - photons_fit_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_trues_frame(blob_index) = positionX_nm_true;
        positionY_nm_trues_frame(blob_index) = positionY_nm_true;
        angleInclination_trues_frame(blob_index) = angleInclination_true;
        angleAzimuth_trues_frame(blob_index) = angleAzimuth_true;
        photons_trues_frame(blob_index) = photons_true;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_fit_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;




        
        
        % % This will create and save a plot for each patch showing the true and estimated positions
        % 
        % % Create figure
        % figure('Position', [100, 100, 600, 500]);
        % 
        % % Display the patch image
        % imagesc(patch_psf_image);
        % colormap('hot'); % Use a colormap that makes the PSF visible
        % axis equal tight;
        % title(sprintf('Frame %d, Blob %d', frame_index, blob_index));
        % colorbar;
        % 
        % % Convert the estimated position from local nm to local pixel coordinates
        % % Note: patch is centered at (half_width+1, half_width+1) in local coordinates
        % local_est_x_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(1) / pixel_size_nm;
        % local_est_y_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(2) / pixel_size_nm;
        % 
        % % Ground truth should be at the center of the patch (approximately)
        % local_gt_x_px = half_width + 1;  % Center of patch
        % local_gt_y_px = half_width + 1;  % Center of patch
        % 
        % hold on;
        % 
        % % Plot true position with a green cross
        % plot(local_gt_x_px, local_gt_y_px, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        % 
        % % Plot estimated position with a red circle
        % plot(local_est_x_px, local_est_y_px, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        % 
        % % Add legend
        % legend('Ground Truth', 'Estimate', 'Location', 'Best');
        % 
        % % Save the figure to the run directory
        % output_filename = sprintf('psf_patch_f%d_b%d.png', frame_index, blob_index);
        % saveas(gcf, fullfile(stack_dir, output_filename));
        % fprintf('Saved visualization to: %s\n', fullfile(stack_dir, output_filename));
        % 
        % % Close the figure
        % close;






        % % DEBUGGING
        % 
        % fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        % fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_true, positionY_nm_true);
        % 
        % % Convert ground truth to pixels
        % gt_x_px = nm_to_px(positionX_nm_true, pixel_size_nm, image_width_px, 'x');
        % gt_y_px = nm_to_px(positionY_nm_true, pixel_size_nm, image_height_px, 'y');
        % fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);
        % 
        % % Show patch info
        % fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
        %         patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        % fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);
        % 
        % % Show actual patch center
        % fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        % fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);
        % 
        % % Calculate ground truth relative to patch center (this should be the target for fitting)
        % fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_true - actual_patch_center_x_nm, ...
        %         positionY_nm_true - actual_patch_center_y_nm);
        % 
        % % After fitting, show the results
        % fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
        %         fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        % fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_estimate, positionY_nm_estimate);
        % fprintf('Error (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_error, positionY_nm_error);
        % fprintf('Error (px): (%.2f, %.2f)\n', ...
        %         positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_trues{frame_index} = positionX_nm_trues_frame;
    positionY_nm_trues{frame_index} = positionY_nm_trues_frame;
    angleInclination_trues{frame_index} = angleInclination_trues_frame;
    angleAzimuth_trues{frame_index} = angleAzimuth_trues_frame;
    photons_trues{frame_index} = photons_trues_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_trues_flat = [positionX_nm_trues{:}];
    positionY_nm_trues_flat = [positionY_nm_trues{:}];
    angleInclination_trues_flat = [angleInclination_trues{:}];
    angleAzimuth_trues_flat = [angleAzimuth_trues{:}];
    photons_trues_flat = [photons_trues{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    photons_errors_flat = [photons_errors{:}];

    % objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_trues_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_trues_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_trues_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_trues_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_trues_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames








%% ----------------------------------------------------------------------
%% HoH
%% ----------------------------------------------------------------------









close all;
clear all;


%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_hinterer_all/";
stack_path = stack_dir + "sim_stack.tif";

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/fitting_results_hinterer_on_hinterer_all.py';

% Model to use
model = 'hinterer';

% Hard-coded microscope parameters
pixel_size_nm = 52;
wavelength = 500;
objectiveNA = 2.17;
objectiveFocalLength = 770;
refractiveIndices = [1.31 , 2.17 , 2.17];
backgroundNoise = 0;
nPhotons = 2000;
pixelSensitivityMask = PixelSensitivity.uniform(9);
nDiscretizationBFP = 129;

patch_width_nm = 988; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
image_width_px = stack_info(1).Width;
image_height_px = stack_info(1).Height;
number_of_frames = numel(stack_info);

% Function to create the current PSF (same as `createFitPSF` method in the class)
% function currentPSF = createFitPSF(psfEstimate, params)
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % Simulate the PSF
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % Apply phase mask
%     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
%     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% 
%     % Apply attenuation mask
%     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
%     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         % Apply aberrations
%         aberrations = getAberrations(psfEstimate,k);
%         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% 
%         % Get image from BFP field
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
%     end
% 
%     currentPsf = adjustExcitation(psfEstimate, currentPsf);
%     currentPsf = applyShotNoise(psfEstimate, currentPsf);
%     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% end
% function currentPSF = createFitPSF(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end
% 
% function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % bfp = BackFocalPlane(psfEstimate);
%     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end




% % Locate all the images and settings files in the dir
% files = dir(frames_dir);
% valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
% valid_extensions_settings = {'.m'};
% frame_paths = {};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_image)
%         frame_paths{end+1} = fullfile(frames_dir, files(i).name);
%     elseif ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(frames_dir, files(i).name);
%     end
% end

% Locate all the settings files in the dir
% files = dir(stack_dir);
% valid_extensions_settings = {'.m'};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(stack_dir, files(i).name);
%     end
% end

files = dir(fullfile(stack_dir, '*.m')); % Use pattern matching
settings_paths = {};
for i = 1:length(files)
    settings_paths{end+1} = fullfile(stack_dir, files(i).name);
end

% Sort the paths based on frame numbers
frame_numbers = zeros(length(settings_paths), 1);
for i = 1:length(settings_paths)
    [~, name, ~] = fileparts(settings_paths{i});
    % Extract the number after 'params_frame'
    frame_number = str2double(regexp(name, 'params_frame(\d+)', 'tokens', 'once'));
    % frame_number = str2double(regexp(name, 'params_(\d+)', 'tokens', 'once'));
    frame_numbers(i) = frame_number;
end

% Sort paths based on frame numbers
[~, idx] = sort(frame_numbers);
settings_paths = settings_paths(idx);

% % Read in thunderstorm localisations
% thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');
% 
% % Extract the required columns
% thunderstorm_frame_array = thunderstorm_results.("frame");
% thunderstorm_x_array = thunderstorm_results.("x [nm]");
% thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % Initialise overall results arrays
% positionX_trues = [];
% positionY_trues = [];
% angleInclination_trues = [];
% angleAzimuth_trues = [];
% 
% positionX_estimates = [];
% positionY_estimates = [];
% angleInclination_estimates = [];
% angleAzimuth_estimates = [];
% 
% positionX_errors = [];
% positionY_errors = [];
% angleInclination_errors = [];
% angleAzimuth_errors = [];
% 
% % objective_function_trues = [];
% % objective_function_estimates = [];
% % objective_function_errors = [];

% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in parameters/settings if using simulated to get ground truth
    evalc('run(settings_path)');
    % par.pixelSize = Length(pixel_size_nm,'nm');
    % par.objectiveNA = objectiveNA;
    % par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    % par.nPixels = min(image_width_px, image_height_px);
    % par.wavelength = Length(wavelength,'nm');
    % 
    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    % par.nDiscretizationBFP = 129;
    % patch_width_nm = 1000; % size of patch around blob to consider
    % patch_width_px = patch_width_nm/pixel_size_nm;
    % patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    % par.nPixels = min(image_width_px, image_height_px);
    % par.nPixels = 2*floor(min(image_width_px, image_height_px)/2) + 1; % This enforces odd number of pixels, as required by Hinterer
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Take thunderstorm results for that frame
    % % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    % current_frame_mask = thunderstorm_frame_array == frame_index;
    % current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    % current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(positionX_nm_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(positionX_nm_array));

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        % adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        % adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);
        % 
        % % Calculate patch boundaries using the adjusted center
        % patch_start_x_px = adjusted_centre_x_px - half_width;
        % patch_start_y_px = adjusted_centre_y_px - half_width;
        % patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        % patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        % 
        % % Extract patch from image
        % patch_indices_x = patch_start_x_px:patch_end_x_px;
        % patch_indices_y = patch_start_y_px:patch_end_y_px;
        % patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        % 
        % % Get actual patch dimensions after boundary checking
        % actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        % actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        % 
        % % Calculate ACTUAL center of patch in pixel coordinates
        % actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        % actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        % 
        % % Convert ACTUAL patch center to physical coordinates
        % actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        % actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        % 
        % % Calculate expected position relative to patch center - beause of
        % % the boundary checking stuff, it might not quite be zero anymore
        % 
        % expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        % expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;





        % Better approach - avoid rounding and maintain sub-pixel precision
        % First calculate the integer patch boundaries that ensure the patch stays within image
        patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
        patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
        patch_end_x_px = min(image_width, patch_start_x_px + patch_width_px - 1);
        patch_end_y_px = min(image_height, patch_start_y_px + patch_width_px - 1);
        
        % Adjust patch start if necessary to maintain constant patch size
        if (patch_end_x_px - patch_start_x_px + 1) < patch_width_px
            patch_start_x_px = max(1, patch_end_x_px - patch_width_px + 1);
        end
        if (patch_end_y_px - patch_start_y_px + 1) < patch_width_px
            patch_start_y_px = max(1, patch_end_y_px - patch_width_px + 1);
        end
        
        % Extract patch from image
        patch_indices_x = patch_start_x_px:patch_end_x_px;
        patch_indices_y = patch_start_y_px:patch_end_y_px;
        patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        
        % Get actual patch dimensions after boundary checking
        actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        
        % Calculate ACTUAL center of patch in pixel coordinates
        actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        
        % Convert ACTUAL patch center to physical coordinates
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        
        % Calculate expected position relative to patch center
        expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;
        





        % % DEBUG output to see coordinates
        % fprintf('Initial patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     patch_centre_x_nm, patch_centre_y_nm, patch_centre_x_px, patch_centre_y_px);
        % fprintf('Actual patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     actual_patch_center_x_nm, actual_patch_center_y_nm, actual_patch_center_x_px, actual_patch_center_y_px);

        % % Clip values just if want to display
        % display_image = (patch_psf_image - min(patch_psf_image(:))) / (max(patch_psf_image(:)) - min(patch_psf_image(:)));
        % imshow(display_image)

        % Setup PSF object for fitting
        par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        if strcmpi(model, 'hinterer')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
            psfInit = PSF(par);
        elseif strcmpi(model, 'mortensen')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_mortensen(par);
        elseif strcmpi(model, 'gaussian')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_gaussian(par);
        else
            error('Unknown model type: %s', model);
        end
        
        psfInit.image = patch_psf_image;
        psfInit.nPixels = length(patch_psf_image);

        % % If not masking or cropping, use this
        % % Generate throwaway PSF object, replace the image with our image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % psfInit.image = psf_image;

        % Run the fit - given nothing, fmincon - reparameterised

        % Give it a background noise estimate (let's say you know it to within 10%)
        if psfInit.backgroundNoise > 1e-5
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly. Will need to change this becaue won't know this for real data
        else
            parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
        end

        parEst.parameterBounds.x = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.newangle1 = [-1, 1];%[1/sqrt(2) - 0.0000000001, 1/sqrt(2) + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];%[0 - 0.0000000001, 0 + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give patch centre as initial value. Add some fake noise to
        % simulate imperfect feature detection.
        parEst.parameterStartValues.x = Length(expected_offset_x + (0 + 200 * (rand() - 0.5)), 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y + (0 + 200 * (rand() - 0.5)), 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
        objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of things like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_true = positionX_nm_array(blob_index);
        positionY_nm_true = positionY_nm_array(blob_index);
        angleInclination_true = angleInclination_array(blob_index);
        angleAzimuth_true = angleAzimuth_array(blob_index);
        photons_true = photon_estimate;

        positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
        angleInclination_error = angleInclination_true - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
        photons_error = photons_true - photons_fit_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_trues_frame(blob_index) = positionX_nm_true;
        positionY_nm_trues_frame(blob_index) = positionY_nm_true;
        angleInclination_trues_frame(blob_index) = angleInclination_true;
        angleAzimuth_trues_frame(blob_index) = angleAzimuth_true;
        photons_trues_frame(blob_index) = photons_true;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_fit_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;




        
        
        % % This will create and save a plot for each patch showing the true and estimated positions
        % 
        % % Create figure
        % figure('Position', [100, 100, 600, 500]);
        % 
        % % Display the patch image
        % imagesc(patch_psf_image);
        % colormap('hot'); % Use a colormap that makes the PSF visible
        % axis equal tight;
        % title(sprintf('Frame %d, Blob %d', frame_index, blob_index));
        % colorbar;
        % 
        % % Convert the estimated position from local nm to local pixel coordinates
        % % Note: patch is centered at (half_width+1, half_width+1) in local coordinates
        % local_est_x_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(1) / pixel_size_nm;
        % local_est_y_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(2) / pixel_size_nm;
        % 
        % % Ground truth should be at the center of the patch (approximately)
        % local_gt_x_px = half_width + 1;  % Center of patch
        % local_gt_y_px = half_width + 1;  % Center of patch
        % 
        % hold on;
        % 
        % % Plot true position with a green cross
        % plot(local_gt_x_px, local_gt_y_px, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        % 
        % % Plot estimated position with a red circle
        % plot(local_est_x_px, local_est_y_px, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        % 
        % % Add legend
        % legend('Ground Truth', 'Estimate', 'Location', 'Best');
        % 
        % % Save the figure to the run directory
        % output_filename = sprintf('psf_patch_f%d_b%d.png', frame_index, blob_index);
        % saveas(gcf, fullfile(stack_dir, output_filename));
        % fprintf('Saved visualization to: %s\n', fullfile(stack_dir, output_filename));
        % 
        % % Close the figure
        % close;






        % % DEBUGGING
        % 
        % fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        % fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_true, positionY_nm_true);
        % 
        % % Convert ground truth to pixels
        % gt_x_px = nm_to_px(positionX_nm_true, pixel_size_nm, image_width_px, 'x');
        % gt_y_px = nm_to_px(positionY_nm_true, pixel_size_nm, image_height_px, 'y');
        % fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);
        % 
        % % Show patch info
        % fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
        %         patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        % fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);
        % 
        % % Show actual patch center
        % fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        % fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);
        % 
        % % Calculate ground truth relative to patch center (this should be the target for fitting)
        % fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_true - actual_patch_center_x_nm, ...
        %         positionY_nm_true - actual_patch_center_y_nm);
        % 
        % % After fitting, show the results
        % fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
        %         fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        % fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_estimate, positionY_nm_estimate);
        % fprintf('Error (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_error, positionY_nm_error);
        % fprintf('Error (px): (%.2f, %.2f)\n', ...
        %         positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_trues{frame_index} = positionX_nm_trues_frame;
    positionY_nm_trues{frame_index} = positionY_nm_trues_frame;
    angleInclination_trues{frame_index} = angleInclination_trues_frame;
    angleAzimuth_trues{frame_index} = angleAzimuth_trues_frame;
    photons_trues{frame_index} = photons_trues_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_trues_flat = [positionX_nm_trues{:}];
    positionY_nm_trues_flat = [positionY_nm_trues{:}];
    angleInclination_trues_flat = [angleInclination_trues{:}];
    angleAzimuth_trues_flat = [angleAzimuth_trues{:}];
    photons_trues_flat = [photons_trues{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    photons_errors_flat = [photons_errors{:}];

    % objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_trues_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_trues_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_trues_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_trues_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_trues_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames








%% ----------------------------------------------------------------------
%% MoH
%% ----------------------------------------------------------------------









close all;
clear all;


%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_hinterer_all/";
stack_path = stack_dir + "sim_stack.tif";

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/fitting_results_mortensen_on_hinterer_all.py';

% Model to use
model = 'mortensen';

% Hard-coded microscope parameters
pixel_size_nm = 52;
wavelength = 500;
objectiveNA = 2.17;
objectiveFocalLength = 770;
refractiveIndices = [1.31 , 2.17 , 2.17];
backgroundNoise = 0;
nPhotons = 2000;
pixelSensitivityMask = PixelSensitivity.uniform(9);
nDiscretizationBFP = 129;

patch_width_nm = 988; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
image_width_px = stack_info(1).Width;
image_height_px = stack_info(1).Height;
number_of_frames = numel(stack_info);

% Function to create the current PSF (same as `createFitPSF` method in the class)
% function currentPSF = createFitPSF(psfEstimate, params)
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % Simulate the PSF
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % Apply phase mask
%     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
%     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% 
%     % Apply attenuation mask
%     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
%     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         % Apply aberrations
%         aberrations = getAberrations(psfEstimate,k);
%         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% 
%         % Get image from BFP field
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
%     end
% 
%     currentPsf = adjustExcitation(psfEstimate, currentPsf);
%     currentPsf = applyShotNoise(psfEstimate, currentPsf);
%     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% end
% function currentPSF = createFitPSF(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end
% 
% function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % bfp = BackFocalPlane(psfEstimate);
%     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end




% % Locate all the images and settings files in the dir
% files = dir(frames_dir);
% valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
% valid_extensions_settings = {'.m'};
% frame_paths = {};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_image)
%         frame_paths{end+1} = fullfile(frames_dir, files(i).name);
%     elseif ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(frames_dir, files(i).name);
%     end
% end

% Locate all the settings files in the dir
% files = dir(stack_dir);
% valid_extensions_settings = {'.m'};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(stack_dir, files(i).name);
%     end
% end

files = dir(fullfile(stack_dir, '*.m')); % Use pattern matching
settings_paths = {};
for i = 1:length(files)
    settings_paths{end+1} = fullfile(stack_dir, files(i).name);
end

% Sort the paths based on frame numbers
frame_numbers = zeros(length(settings_paths), 1);
for i = 1:length(settings_paths)
    [~, name, ~] = fileparts(settings_paths{i});
    % Extract the number after 'params_frame'
    frame_number = str2double(regexp(name, 'params_frame(\d+)', 'tokens', 'once'));
    % frame_number = str2double(regexp(name, 'params_(\d+)', 'tokens', 'once'));
    frame_numbers(i) = frame_number;
end

% Sort paths based on frame numbers
[~, idx] = sort(frame_numbers);
settings_paths = settings_paths(idx);

% % Read in thunderstorm localisations
% thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');
% 
% % Extract the required columns
% thunderstorm_frame_array = thunderstorm_results.("frame");
% thunderstorm_x_array = thunderstorm_results.("x [nm]");
% thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % Initialise overall results arrays
% positionX_trues = [];
% positionY_trues = [];
% angleInclination_trues = [];
% angleAzimuth_trues = [];
% 
% positionX_estimates = [];
% positionY_estimates = [];
% angleInclination_estimates = [];
% angleAzimuth_estimates = [];
% 
% positionX_errors = [];
% positionY_errors = [];
% angleInclination_errors = [];
% angleAzimuth_errors = [];
% 
% % objective_function_trues = [];
% % objective_function_estimates = [];
% % objective_function_errors = [];

% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in parameters/settings if using simulated to get ground truth
    evalc('run(settings_path)');
    % par.pixelSize = Length(pixel_size_nm,'nm');
    % par.objectiveNA = objectiveNA;
    % par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    % par.nPixels = min(image_width_px, image_height_px);
    % par.wavelength = Length(wavelength,'nm');
    % 
    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    % par.nDiscretizationBFP = 129;
    % patch_width_nm = 1000; % size of patch around blob to consider
    % patch_width_px = patch_width_nm/pixel_size_nm;
    % patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    % par.nPixels = min(image_width_px, image_height_px);
    % par.nPixels = 2*floor(min(image_width_px, image_height_px)/2) + 1; % This enforces odd number of pixels, as required by Hinterer
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Take thunderstorm results for that frame
    % % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    % current_frame_mask = thunderstorm_frame_array == frame_index;
    % current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    % current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(positionX_nm_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(positionX_nm_array));

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        % adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        % adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);
        % 
        % % Calculate patch boundaries using the adjusted center
        % patch_start_x_px = adjusted_centre_x_px - half_width;
        % patch_start_y_px = adjusted_centre_y_px - half_width;
        % patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        % patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        % 
        % % Extract patch from image
        % patch_indices_x = patch_start_x_px:patch_end_x_px;
        % patch_indices_y = patch_start_y_px:patch_end_y_px;
        % patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        % 
        % % Get actual patch dimensions after boundary checking
        % actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        % actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        % 
        % % Calculate ACTUAL center of patch in pixel coordinates
        % actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        % actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        % 
        % % Convert ACTUAL patch center to physical coordinates
        % actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        % actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        % 
        % % Calculate expected position relative to patch center - beause of
        % % the boundary checking stuff, it might not quite be zero anymore
        % 
        % expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        % expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;





        % Better approach - avoid rounding and maintain sub-pixel precision
        % First calculate the integer patch boundaries that ensure the patch stays within image
        patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
        patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
        patch_end_x_px = min(image_width, patch_start_x_px + patch_width_px - 1);
        patch_end_y_px = min(image_height, patch_start_y_px + patch_width_px - 1);
        
        % Adjust patch start if necessary to maintain constant patch size
        if (patch_end_x_px - patch_start_x_px + 1) < patch_width_px
            patch_start_x_px = max(1, patch_end_x_px - patch_width_px + 1);
        end
        if (patch_end_y_px - patch_start_y_px + 1) < patch_width_px
            patch_start_y_px = max(1, patch_end_y_px - patch_width_px + 1);
        end
        
        % Extract patch from image
        patch_indices_x = patch_start_x_px:patch_end_x_px;
        patch_indices_y = patch_start_y_px:patch_end_y_px;
        patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        
        % Get actual patch dimensions after boundary checking
        actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        
        % Calculate ACTUAL center of patch in pixel coordinates
        actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        
        % Convert ACTUAL patch center to physical coordinates
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        
        % Calculate expected position relative to patch center
        expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;
        





        % % DEBUG output to see coordinates
        % fprintf('Initial patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     patch_centre_x_nm, patch_centre_y_nm, patch_centre_x_px, patch_centre_y_px);
        % fprintf('Actual patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     actual_patch_center_x_nm, actual_patch_center_y_nm, actual_patch_center_x_px, actual_patch_center_y_px);

        % % Clip values just if want to display
        % display_image = (patch_psf_image - min(patch_psf_image(:))) / (max(patch_psf_image(:)) - min(patch_psf_image(:)));
        % imshow(display_image)

        % Setup PSF object for fitting
        par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        if strcmpi(model, 'hinterer')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
            psfInit = PSF(par);
        elseif strcmpi(model, 'mortensen')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_mortensen(par);
        elseif strcmpi(model, 'gaussian')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_gaussian(par);
        else
            error('Unknown model type: %s', model);
        end
        
        psfInit.image = patch_psf_image;
        psfInit.nPixels = length(patch_psf_image);

        % % If not masking or cropping, use this
        % % Generate throwaway PSF object, replace the image with our image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % psfInit.image = psf_image;

        % Run the fit - given nothing, fmincon - reparameterised

        % Give it a background noise estimate (let's say you know it to within 10%)
        if psfInit.backgroundNoise > 1e-5
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly. Will need to change this becaue won't know this for real data
        else
            parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
        end

        parEst.parameterBounds.x = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.newangle1 = [-1, 1];%[1/sqrt(2) - 0.0000000001, 1/sqrt(2) + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];%[0 - 0.0000000001, 0 + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give patch centre as initial value. Add some fake noise to
        % simulate imperfect feature detection.
        parEst.parameterStartValues.x = Length(expected_offset_x + (0 + 200 * (rand() - 0.5)), 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y + (0 + 200 * (rand() - 0.5)), 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
        objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of things like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_true = positionX_nm_array(blob_index);
        positionY_nm_true = positionY_nm_array(blob_index);
        angleInclination_true = angleInclination_array(blob_index);
        angleAzimuth_true = angleAzimuth_array(blob_index);
        photons_true = photon_estimate;

        positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
        angleInclination_error = angleInclination_true - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
        photons_error = photons_true - photons_fit_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_trues_frame(blob_index) = positionX_nm_true;
        positionY_nm_trues_frame(blob_index) = positionY_nm_true;
        angleInclination_trues_frame(blob_index) = angleInclination_true;
        angleAzimuth_trues_frame(blob_index) = angleAzimuth_true;
        photons_trues_frame(blob_index) = photons_true;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_fit_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;




        
        
        % % This will create and save a plot for each patch showing the true and estimated positions
        % 
        % % Create figure
        % figure('Position', [100, 100, 600, 500]);
        % 
        % % Display the patch image
        % imagesc(patch_psf_image);
        % colormap('hot'); % Use a colormap that makes the PSF visible
        % axis equal tight;
        % title(sprintf('Frame %d, Blob %d', frame_index, blob_index));
        % colorbar;
        % 
        % % Convert the estimated position from local nm to local pixel coordinates
        % % Note: patch is centered at (half_width+1, half_width+1) in local coordinates
        % local_est_x_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(1) / pixel_size_nm;
        % local_est_y_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(2) / pixel_size_nm;
        % 
        % % Ground truth should be at the center of the patch (approximately)
        % local_gt_x_px = half_width + 1;  % Center of patch
        % local_gt_y_px = half_width + 1;  % Center of patch
        % 
        % hold on;
        % 
        % % Plot true position with a green cross
        % plot(local_gt_x_px, local_gt_y_px, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        % 
        % % Plot estimated position with a red circle
        % plot(local_est_x_px, local_est_y_px, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        % 
        % % Add legend
        % legend('Ground Truth', 'Estimate', 'Location', 'Best');
        % 
        % % Save the figure to the run directory
        % output_filename = sprintf('psf_patch_f%d_b%d.png', frame_index, blob_index);
        % saveas(gcf, fullfile(stack_dir, output_filename));
        % fprintf('Saved visualization to: %s\n', fullfile(stack_dir, output_filename));
        % 
        % % Close the figure
        % close;






        % % DEBUGGING
        % 
        % fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        % fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_true, positionY_nm_true);
        % 
        % % Convert ground truth to pixels
        % gt_x_px = nm_to_px(positionX_nm_true, pixel_size_nm, image_width_px, 'x');
        % gt_y_px = nm_to_px(positionY_nm_true, pixel_size_nm, image_height_px, 'y');
        % fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);
        % 
        % % Show patch info
        % fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
        %         patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        % fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);
        % 
        % % Show actual patch center
        % fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        % fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);
        % 
        % % Calculate ground truth relative to patch center (this should be the target for fitting)
        % fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_true - actual_patch_center_x_nm, ...
        %         positionY_nm_true - actual_patch_center_y_nm);
        % 
        % % After fitting, show the results
        % fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
        %         fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        % fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_estimate, positionY_nm_estimate);
        % fprintf('Error (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_error, positionY_nm_error);
        % fprintf('Error (px): (%.2f, %.2f)\n', ...
        %         positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_trues{frame_index} = positionX_nm_trues_frame;
    positionY_nm_trues{frame_index} = positionY_nm_trues_frame;
    angleInclination_trues{frame_index} = angleInclination_trues_frame;
    angleAzimuth_trues{frame_index} = angleAzimuth_trues_frame;
    photons_trues{frame_index} = photons_trues_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_trues_flat = [positionX_nm_trues{:}];
    positionY_nm_trues_flat = [positionY_nm_trues{:}];
    angleInclination_trues_flat = [angleInclination_trues{:}];
    angleAzimuth_trues_flat = [angleAzimuth_trues{:}];
    photons_trues_flat = [photons_trues{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    photons_errors_flat = [photons_errors{:}];

    % objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_trues_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_trues_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_trues_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_trues_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_trues_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames








%% ----------------------------------------------------------------------
%% GoM
%% ----------------------------------------------------------------------









close all;
clear all;


%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_all/";
stack_path = stack_dir + "sim_stack.tif";

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/fitting_results_gaussian_on_mortensen_all.py';

% Model to use
model = 'gaussian';

% Hard-coded microscope parameters
pixel_size_nm = 52;
wavelength = 500;
objectiveNA = 2.17;
objectiveFocalLength = 770;
refractiveIndices = [1.31 , 2.17 , 2.17];
backgroundNoise = 0;
nPhotons = 2000;
pixelSensitivityMask = PixelSensitivity.uniform(9);
nDiscretizationBFP = 129;

patch_width_nm = 988; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
image_width_px = stack_info(1).Width;
image_height_px = stack_info(1).Height;
number_of_frames = numel(stack_info);

% Function to create the current PSF (same as `createFitPSF` method in the class)
% function currentPSF = createFitPSF(psfEstimate, params)
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % Simulate the PSF
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % Apply phase mask
%     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
%     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% 
%     % Apply attenuation mask
%     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
%     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         % Apply aberrations
%         aberrations = getAberrations(psfEstimate,k);
%         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% 
%         % Get image from BFP field
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
%     end
% 
%     currentPsf = adjustExcitation(psfEstimate, currentPsf);
%     currentPsf = applyShotNoise(psfEstimate, currentPsf);
%     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% end
% function currentPSF = createFitPSF(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end
% 
% function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % bfp = BackFocalPlane(psfEstimate);
%     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end




% % Locate all the images and settings files in the dir
% files = dir(frames_dir);
% valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
% valid_extensions_settings = {'.m'};
% frame_paths = {};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_image)
%         frame_paths{end+1} = fullfile(frames_dir, files(i).name);
%     elseif ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(frames_dir, files(i).name);
%     end
% end

% Locate all the settings files in the dir
% files = dir(stack_dir);
% valid_extensions_settings = {'.m'};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(stack_dir, files(i).name);
%     end
% end

files = dir(fullfile(stack_dir, '*.m')); % Use pattern matching
settings_paths = {};
for i = 1:length(files)
    settings_paths{end+1} = fullfile(stack_dir, files(i).name);
end

% Sort the paths based on frame numbers
frame_numbers = zeros(length(settings_paths), 1);
for i = 1:length(settings_paths)
    [~, name, ~] = fileparts(settings_paths{i});
    % Extract the number after 'params_frame'
    frame_number = str2double(regexp(name, 'params_frame(\d+)', 'tokens', 'once'));
    % frame_number = str2double(regexp(name, 'params_(\d+)', 'tokens', 'once'));
    frame_numbers(i) = frame_number;
end

% Sort paths based on frame numbers
[~, idx] = sort(frame_numbers);
settings_paths = settings_paths(idx);

% % Read in thunderstorm localisations
% thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');
% 
% % Extract the required columns
% thunderstorm_frame_array = thunderstorm_results.("frame");
% thunderstorm_x_array = thunderstorm_results.("x [nm]");
% thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % Initialise overall results arrays
% positionX_trues = [];
% positionY_trues = [];
% angleInclination_trues = [];
% angleAzimuth_trues = [];
% 
% positionX_estimates = [];
% positionY_estimates = [];
% angleInclination_estimates = [];
% angleAzimuth_estimates = [];
% 
% positionX_errors = [];
% positionY_errors = [];
% angleInclination_errors = [];
% angleAzimuth_errors = [];
% 
% % objective_function_trues = [];
% % objective_function_estimates = [];
% % objective_function_errors = [];

% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in parameters/settings if using simulated to get ground truth
    evalc('run(settings_path)');
    % par.pixelSize = Length(pixel_size_nm,'nm');
    % par.objectiveNA = objectiveNA;
    % par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    % par.nPixels = min(image_width_px, image_height_px);
    % par.wavelength = Length(wavelength,'nm');
    % 
    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    % par.nDiscretizationBFP = 129;
    % patch_width_nm = 1000; % size of patch around blob to consider
    % patch_width_px = patch_width_nm/pixel_size_nm;
    % patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    % par.nPixels = min(image_width_px, image_height_px);
    % par.nPixels = 2*floor(min(image_width_px, image_height_px)/2) + 1; % This enforces odd number of pixels, as required by Hinterer
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Take thunderstorm results for that frame
    % % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    % current_frame_mask = thunderstorm_frame_array == frame_index;
    % current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    % current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(positionX_nm_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(positionX_nm_array));

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        % adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        % adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);
        % 
        % % Calculate patch boundaries using the adjusted center
        % patch_start_x_px = adjusted_centre_x_px - half_width;
        % patch_start_y_px = adjusted_centre_y_px - half_width;
        % patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        % patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        % 
        % % Extract patch from image
        % patch_indices_x = patch_start_x_px:patch_end_x_px;
        % patch_indices_y = patch_start_y_px:patch_end_y_px;
        % patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        % 
        % % Get actual patch dimensions after boundary checking
        % actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        % actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        % 
        % % Calculate ACTUAL center of patch in pixel coordinates
        % actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        % actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        % 
        % % Convert ACTUAL patch center to physical coordinates
        % actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        % actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        % 
        % % Calculate expected position relative to patch center - beause of
        % % the boundary checking stuff, it might not quite be zero anymore
        % 
        % expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        % expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;





        % Better approach - avoid rounding and maintain sub-pixel precision
        % First calculate the integer patch boundaries that ensure the patch stays within image
        patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
        patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
        patch_end_x_px = min(image_width, patch_start_x_px + patch_width_px - 1);
        patch_end_y_px = min(image_height, patch_start_y_px + patch_width_px - 1);
        
        % Adjust patch start if necessary to maintain constant patch size
        if (patch_end_x_px - patch_start_x_px + 1) < patch_width_px
            patch_start_x_px = max(1, patch_end_x_px - patch_width_px + 1);
        end
        if (patch_end_y_px - patch_start_y_px + 1) < patch_width_px
            patch_start_y_px = max(1, patch_end_y_px - patch_width_px + 1);
        end
        
        % Extract patch from image
        patch_indices_x = patch_start_x_px:patch_end_x_px;
        patch_indices_y = patch_start_y_px:patch_end_y_px;
        patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        
        % Get actual patch dimensions after boundary checking
        actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        
        % Calculate ACTUAL center of patch in pixel coordinates
        actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        
        % Convert ACTUAL patch center to physical coordinates
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        
        % Calculate expected position relative to patch center
        expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;
        





        % % DEBUG output to see coordinates
        % fprintf('Initial patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     patch_centre_x_nm, patch_centre_y_nm, patch_centre_x_px, patch_centre_y_px);
        % fprintf('Actual patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     actual_patch_center_x_nm, actual_patch_center_y_nm, actual_patch_center_x_px, actual_patch_center_y_px);

        % % Clip values just if want to display
        % display_image = (patch_psf_image - min(patch_psf_image(:))) / (max(patch_psf_image(:)) - min(patch_psf_image(:)));
        % imshow(display_image)

        % Setup PSF object for fitting
        par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        if strcmpi(model, 'hinterer')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
            psfInit = PSF(par);
        elseif strcmpi(model, 'mortensen')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_mortensen(par);
        elseif strcmpi(model, 'gaussian')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_gaussian(par);
        else
            error('Unknown model type: %s', model);
        end
        
        psfInit.image = patch_psf_image;
        psfInit.nPixels = length(patch_psf_image);

        % % If not masking or cropping, use this
        % % Generate throwaway PSF object, replace the image with our image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % psfInit.image = psf_image;

        % Run the fit - given nothing, fmincon - reparameterised

        % Give it a background noise estimate (let's say you know it to within 10%)
        if psfInit.backgroundNoise > 1e-5
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly. Will need to change this becaue won't know this for real data
        else
            parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
        end

        parEst.parameterBounds.x = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.newangle1 = [-1, 1];%[1/sqrt(2) - 0.0000000001, 1/sqrt(2) + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];%[0 - 0.0000000001, 0 + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give patch centre as initial value. Add some fake noise to
        % simulate imperfect feature detection.
        parEst.parameterStartValues.x = Length(expected_offset_x + (0 + 200 * (rand() - 0.5)), 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y + (0 + 200 * (rand() - 0.5)), 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
        objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of things like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_true = positionX_nm_array(blob_index);
        positionY_nm_true = positionY_nm_array(blob_index);
        angleInclination_true = angleInclination_array(blob_index);
        angleAzimuth_true = angleAzimuth_array(blob_index);
        photons_true = photon_estimate;

        positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
        angleInclination_error = angleInclination_true - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
        photons_error = photons_true - photons_fit_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_trues_frame(blob_index) = positionX_nm_true;
        positionY_nm_trues_frame(blob_index) = positionY_nm_true;
        angleInclination_trues_frame(blob_index) = angleInclination_true;
        angleAzimuth_trues_frame(blob_index) = angleAzimuth_true;
        photons_trues_frame(blob_index) = photons_true;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_fit_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;




        
        
        % % This will create and save a plot for each patch showing the true and estimated positions
        % 
        % % Create figure
        % figure('Position', [100, 100, 600, 500]);
        % 
        % % Display the patch image
        % imagesc(patch_psf_image);
        % colormap('hot'); % Use a colormap that makes the PSF visible
        % axis equal tight;
        % title(sprintf('Frame %d, Blob %d', frame_index, blob_index));
        % colorbar;
        % 
        % % Convert the estimated position from local nm to local pixel coordinates
        % % Note: patch is centered at (half_width+1, half_width+1) in local coordinates
        % local_est_x_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(1) / pixel_size_nm;
        % local_est_y_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(2) / pixel_size_nm;
        % 
        % % Ground truth should be at the center of the patch (approximately)
        % local_gt_x_px = half_width + 1;  % Center of patch
        % local_gt_y_px = half_width + 1;  % Center of patch
        % 
        % hold on;
        % 
        % % Plot true position with a green cross
        % plot(local_gt_x_px, local_gt_y_px, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        % 
        % % Plot estimated position with a red circle
        % plot(local_est_x_px, local_est_y_px, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        % 
        % % Add legend
        % legend('Ground Truth', 'Estimate', 'Location', 'Best');
        % 
        % % Save the figure to the run directory
        % output_filename = sprintf('psf_patch_f%d_b%d.png', frame_index, blob_index);
        % saveas(gcf, fullfile(stack_dir, output_filename));
        % fprintf('Saved visualization to: %s\n', fullfile(stack_dir, output_filename));
        % 
        % % Close the figure
        % close;






        % % DEBUGGING
        % 
        % fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        % fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_true, positionY_nm_true);
        % 
        % % Convert ground truth to pixels
        % gt_x_px = nm_to_px(positionX_nm_true, pixel_size_nm, image_width_px, 'x');
        % gt_y_px = nm_to_px(positionY_nm_true, pixel_size_nm, image_height_px, 'y');
        % fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);
        % 
        % % Show patch info
        % fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
        %         patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        % fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);
        % 
        % % Show actual patch center
        % fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        % fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);
        % 
        % % Calculate ground truth relative to patch center (this should be the target for fitting)
        % fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_true - actual_patch_center_x_nm, ...
        %         positionY_nm_true - actual_patch_center_y_nm);
        % 
        % % After fitting, show the results
        % fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
        %         fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        % fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_estimate, positionY_nm_estimate);
        % fprintf('Error (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_error, positionY_nm_error);
        % fprintf('Error (px): (%.2f, %.2f)\n', ...
        %         positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_trues{frame_index} = positionX_nm_trues_frame;
    positionY_nm_trues{frame_index} = positionY_nm_trues_frame;
    angleInclination_trues{frame_index} = angleInclination_trues_frame;
    angleAzimuth_trues{frame_index} = angleAzimuth_trues_frame;
    photons_trues{frame_index} = photons_trues_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_trues_flat = [positionX_nm_trues{:}];
    positionY_nm_trues_flat = [positionY_nm_trues{:}];
    angleInclination_trues_flat = [angleInclination_trues{:}];
    angleAzimuth_trues_flat = [angleAzimuth_trues{:}];
    photons_trues_flat = [photons_trues{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    photons_errors_flat = [photons_errors{:}];

    % objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_trues_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_trues_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_trues_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_trues_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_trues_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames








%% ----------------------------------------------------------------------
%% GoH
%% ----------------------------------------------------------------------







close all;
clear all;


%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_hinterer_all/";
stack_path = stack_dir + "sim_stack.tif";

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/fitting_results_gaussian_on_hinterer_all.py';

% Model to use
model = 'gaussian';

% Hard-coded microscope parameters
pixel_size_nm = 52;
wavelength = 500;
objectiveNA = 2.17;
objectiveFocalLength = 770;
refractiveIndices = [1.31 , 2.17 , 2.17];
backgroundNoise = 0;
nPhotons = 2000;
pixelSensitivityMask = PixelSensitivity.uniform(9);
nDiscretizationBFP = 129;

patch_width_nm = 988; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
image_width_px = stack_info(1).Width;
image_height_px = stack_info(1).Height;
number_of_frames = numel(stack_info);

% Function to create the current PSF (same as `createFitPSF` method in the class)
% function currentPSF = createFitPSF(psfEstimate, params)
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % Simulate the PSF
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % Apply phase mask
%     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
%     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% 
%     % Apply attenuation mask
%     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
%     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         % Apply aberrations
%         aberrations = getAberrations(psfEstimate,k);
%         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% 
%         % Get image from BFP field
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
%     end
% 
%     currentPsf = adjustExcitation(psfEstimate, currentPsf);
%     currentPsf = applyShotNoise(psfEstimate, currentPsf);
%     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% end
% function currentPSF = createFitPSF(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     bfp = BackFocalPlane(psfEstimate);
%     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end
% 
% function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% 
%     psfEstimate.position = Length([params(1:2), 0], 'nm');
%     psfEstimate.defocus = Length(params(3), 'nm');
%     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
%     noiseEstimate = psfEstimate.backgroundNoise;
%     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% 
%     % bfp = BackFocalPlane(psfEstimate);
%     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
%     psfEstimate.backFocalPlane = bfp;
% 
%     % dave
%     psfEstimate.fieldBFP.x = bfp.electricField.x;
%     psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%     for k=1:size(psfEstimate.stageDrift.motion,1)
%         aberrationCoeffs = getAberrations(psfEstimate,k);
%         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%     end
% 
%     totalIntensity = sum(currentPsf,'all');
%     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
%     currentFitPSF = currentPsf ./ norm(currentPsf);
% 
%     currentPSF = currentFitPSF;
% 
% end




% % Locate all the images and settings files in the dir
% files = dir(frames_dir);
% valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
% valid_extensions_settings = {'.m'};
% frame_paths = {};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_image)
%         frame_paths{end+1} = fullfile(frames_dir, files(i).name);
%     elseif ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(frames_dir, files(i).name);
%     end
% end

% Locate all the settings files in the dir
% files = dir(stack_dir);
% valid_extensions_settings = {'.m'};
% settings_paths = {};
% 
% for i = 1:length(files)
%     [~, ~, ext] = fileparts(files(i).name);
%     if ismember(lower(ext), valid_extensions_settings)
%         settings_paths{end+1} = fullfile(stack_dir, files(i).name);
%     end
% end

files = dir(fullfile(stack_dir, '*.m')); % Use pattern matching
settings_paths = {};
for i = 1:length(files)
    settings_paths{end+1} = fullfile(stack_dir, files(i).name);
end

% Sort the paths based on frame numbers
frame_numbers = zeros(length(settings_paths), 1);
for i = 1:length(settings_paths)
    [~, name, ~] = fileparts(settings_paths{i});
    % Extract the number after 'params_frame'
    frame_number = str2double(regexp(name, 'params_frame(\d+)', 'tokens', 'once'));
    % frame_number = str2double(regexp(name, 'params_(\d+)', 'tokens', 'once'));
    frame_numbers(i) = frame_number;
end

% Sort paths based on frame numbers
[~, idx] = sort(frame_numbers);
settings_paths = settings_paths(idx);

% % Read in thunderstorm localisations
% thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');
% 
% % Extract the required columns
% thunderstorm_frame_array = thunderstorm_results.("frame");
% thunderstorm_x_array = thunderstorm_results.("x [nm]");
% thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % Initialise overall results arrays
% positionX_trues = [];
% positionY_trues = [];
% angleInclination_trues = [];
% angleAzimuth_trues = [];
% 
% positionX_estimates = [];
% positionY_estimates = [];
% angleInclination_estimates = [];
% angleAzimuth_estimates = [];
% 
% positionX_errors = [];
% positionY_errors = [];
% angleInclination_errors = [];
% angleAzimuth_errors = [];
% 
% % objective_function_trues = [];
% % objective_function_estimates = [];
% % objective_function_errors = [];

% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in parameters/settings if using simulated to get ground truth
    evalc('run(settings_path)');
    % par.pixelSize = Length(pixel_size_nm,'nm');
    % par.objectiveNA = objectiveNA;
    % par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    % par.nPixels = min(image_width_px, image_height_px);
    % par.wavelength = Length(wavelength,'nm');
    % 
    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    % par.nDiscretizationBFP = 129;
    % patch_width_nm = 1000; % size of patch around blob to consider
    % patch_width_px = patch_width_nm/pixel_size_nm;
    % patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    % par.nPixels = min(image_width_px, image_height_px);
    % par.nPixels = 2*floor(min(image_width_px, image_height_px)/2) + 1; % This enforces odd number of pixels, as required by Hinterer
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Take thunderstorm results for that frame
    % % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    % current_frame_mask = thunderstorm_frame_array == frame_index;
    % current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    % current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(positionX_nm_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(positionX_nm_array));

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        % adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        % adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);
        % 
        % % Calculate patch boundaries using the adjusted center
        % patch_start_x_px = adjusted_centre_x_px - half_width;
        % patch_start_y_px = adjusted_centre_y_px - half_width;
        % patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        % patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        % 
        % % Extract patch from image
        % patch_indices_x = patch_start_x_px:patch_end_x_px;
        % patch_indices_y = patch_start_y_px:patch_end_y_px;
        % patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        % 
        % % Get actual patch dimensions after boundary checking
        % actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        % actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        % 
        % % Calculate ACTUAL center of patch in pixel coordinates
        % actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        % actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        % 
        % % Convert ACTUAL patch center to physical coordinates
        % actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        % actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        % 
        % % Calculate expected position relative to patch center - beause of
        % % the boundary checking stuff, it might not quite be zero anymore
        % 
        % expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        % expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;





        % Better approach - avoid rounding and maintain sub-pixel precision
        % First calculate the integer patch boundaries that ensure the patch stays within image
        patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
        patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
        patch_end_x_px = min(image_width, patch_start_x_px + patch_width_px - 1);
        patch_end_y_px = min(image_height, patch_start_y_px + patch_width_px - 1);
        
        % Adjust patch start if necessary to maintain constant patch size
        if (patch_end_x_px - patch_start_x_px + 1) < patch_width_px
            patch_start_x_px = max(1, patch_end_x_px - patch_width_px + 1);
        end
        if (patch_end_y_px - patch_start_y_px + 1) < patch_width_px
            patch_start_y_px = max(1, patch_end_y_px - patch_width_px + 1);
        end
        
        % Extract patch from image
        patch_indices_x = patch_start_x_px:patch_end_x_px;
        patch_indices_y = patch_start_y_px:patch_end_y_px;
        patch_psf_image = psf_image(patch_indices_y, patch_indices_x);
        
        % Get actual patch dimensions after boundary checking
        actual_patch_width_px = patch_end_x_px - patch_start_x_px + 1;
        actual_patch_height_px = patch_end_y_px - patch_start_y_px + 1;
        
        % Calculate ACTUAL center of patch in pixel coordinates
        actual_patch_center_x_px = patch_start_x_px + (actual_patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (actual_patch_height_px - 1)/2;
        
        % Convert ACTUAL patch center to physical coordinates
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_width_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_height_px, 'y');
        
        % Calculate expected position relative to patch center
        expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
        expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;
        





        % % DEBUG output to see coordinates
        % fprintf('Initial patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     patch_centre_x_nm, patch_centre_y_nm, patch_centre_x_px, patch_centre_y_px);
        % fprintf('Actual patch center: (%.2f, %.2f) nm, (%.2f, %.2f) px\n', ...
        %     actual_patch_center_x_nm, actual_patch_center_y_nm, actual_patch_center_x_px, actual_patch_center_y_px);

        % % Clip values just if want to display
        % display_image = (patch_psf_image - min(patch_psf_image(:))) / (max(patch_psf_image(:)) - min(patch_psf_image(:)));
        % imshow(display_image)

        % Setup PSF object for fitting
        par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        if strcmpi(model, 'hinterer')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
            psfInit = PSF(par);
        elseif strcmpi(model, 'mortensen')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_mortensen(par);
        elseif strcmpi(model, 'gaussian')
            % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
            psfInit = PSF_gaussian(par);
        else
            error('Unknown model type: %s', model);
        end
        
        psfInit.image = patch_psf_image;
        psfInit.nPixels = length(patch_psf_image);

        % % If not masking or cropping, use this
        % % Generate throwaway PSF object, replace the image with our image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % psfInit.image = psf_image;

        % Run the fit - given nothing, fmincon - reparameterised

        % Give it a background noise estimate (let's say you know it to within 10%)
        if psfInit.backgroundNoise > 1e-5
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly. Will need to change this becaue won't know this for real data
        else
            parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
        end

        parEst.parameterBounds.x = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-patch_width_nm/2 patch_width_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.newangle1 = [-1, 1];%[1/sqrt(2) - 0.0000000001, 1/sqrt(2) + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];%[0 - 0.0000000001, 0 + 0.0000000001];%[-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give patch centre as initial value. Add some fake noise to
        % simulate imperfect feature detection.
        parEst.parameterStartValues.x = Length(expected_offset_x + (0 + 200 * (rand() - 0.5)), 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y + (0 + 200 * (rand() - 0.5)), 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
        objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of things like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_true = positionX_nm_array(blob_index);
        positionY_nm_true = positionY_nm_array(blob_index);
        angleInclination_true = angleInclination_array(blob_index);
        angleAzimuth_true = angleAzimuth_array(blob_index);
        photons_true = photon_estimate;

        positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
        angleInclination_error = angleInclination_true - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
        photons_error = photons_true - photons_fit_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_trues_frame(blob_index) = positionX_nm_true;
        positionY_nm_trues_frame(blob_index) = positionY_nm_true;
        angleInclination_trues_frame(blob_index) = angleInclination_true;
        angleAzimuth_trues_frame(blob_index) = angleAzimuth_true;
        photons_trues_frame(blob_index) = photons_true;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_fit_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;




        
        
        % % This will create and save a plot for each patch showing the true and estimated positions
        % 
        % % Create figure
        % figure('Position', [100, 100, 600, 500]);
        % 
        % % Display the patch image
        % imagesc(patch_psf_image);
        % colormap('hot'); % Use a colormap that makes the PSF visible
        % axis equal tight;
        % title(sprintf('Frame %d, Blob %d', frame_index, blob_index));
        % colorbar;
        % 
        % % Convert the estimated position from local nm to local pixel coordinates
        % % Note: patch is centered at (half_width+1, half_width+1) in local coordinates
        % local_est_x_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(1) / pixel_size_nm;
        % local_est_y_px = half_width + 1 + fitResult.estimatesPositionDefocus.ML(2) / pixel_size_nm;
        % 
        % % Ground truth should be at the center of the patch (approximately)
        % local_gt_x_px = half_width + 1;  % Center of patch
        % local_gt_y_px = half_width + 1;  % Center of patch
        % 
        % hold on;
        % 
        % % Plot true position with a green cross
        % plot(local_gt_x_px, local_gt_y_px, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        % 
        % % Plot estimated position with a red circle
        % plot(local_est_x_px, local_est_y_px, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        % 
        % % Add legend
        % legend('Ground Truth', 'Estimate', 'Location', 'Best');
        % 
        % % Save the figure to the run directory
        % output_filename = sprintf('psf_patch_f%d_b%d.png', frame_index, blob_index);
        % saveas(gcf, fullfile(stack_dir, output_filename));
        % fprintf('Saved visualization to: %s\n', fullfile(stack_dir, output_filename));
        % 
        % % Close the figure
        % close;






        % % DEBUGGING
        % 
        % fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        % fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_true, positionY_nm_true);
        % 
        % % Convert ground truth to pixels
        % gt_x_px = nm_to_px(positionX_nm_true, pixel_size_nm, image_width_px, 'x');
        % gt_y_px = nm_to_px(positionY_nm_true, pixel_size_nm, image_height_px, 'y');
        % fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);
        % 
        % % Show patch info
        % fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
        %         patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        % fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);
        % 
        % % Show actual patch center
        % fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        % fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);
        % 
        % % Calculate ground truth relative to patch center (this should be the target for fitting)
        % fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_true - actual_patch_center_x_nm, ...
        %         positionY_nm_true - actual_patch_center_y_nm);
        % 
        % % After fitting, show the results
        % fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
        %         fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        % fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_estimate, positionY_nm_estimate);
        % fprintf('Error (nm): (%.2f, %.2f)\n', ...
        %         positionX_nm_error, positionY_nm_error);
        % fprintf('Error (px): (%.2f, %.2f)\n', ...
        %         positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_trues{frame_index} = positionX_nm_trues_frame;
    positionY_nm_trues{frame_index} = positionY_nm_trues_frame;
    angleInclination_trues{frame_index} = angleInclination_trues_frame;
    angleAzimuth_trues{frame_index} = angleAzimuth_trues_frame;
    photons_trues{frame_index} = photons_trues_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_trues_flat = [positionX_nm_trues{:}];
    positionY_nm_trues_flat = [positionY_nm_trues{:}];
    angleInclination_trues_flat = [angleInclination_trues{:}];
    angleAzimuth_trues_flat = [angleAzimuth_trues{:}];
    photons_trues_flat = [photons_trues{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    photons_errors_flat = [photons_errors{:}];

    % objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_trues_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_trues_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_trues_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_trues_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_trues_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames
