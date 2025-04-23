
%% Thunderstorm patch version


%% Fitting multiple PSFs in a single frame

% Same as loop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Fit
%% ----------

% Input params
stack_dir = "/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/real/";
stack_path = stack_dir + "stack_all_simulations.tif";

% Thunderstorm localisations
thunderstorm_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/real/thunderstorm_results.csv';

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/real/results_real.py';

% Model to use
model = 'hinterer';

% ROI location in nm, with origin at centre of image
ROI_centre_x = 0;%1234;
ROI_centre_y = 0;%2*ROI_centre_x;
ROI_width = 1000000;%4000;
ROI_height = ROI_width/2;
ROI_min_x = ROI_centre_x - ROI_width/2;
ROI_max_x = ROI_centre_x + ROI_width/2;
ROI_min_y = ROI_centre_y - ROI_height/2;
ROI_max_y = ROI_centre_y + ROI_height/2;

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

patch_width_nm = 1000; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;
patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

% Pull out info about image stack
stack_info = imfinfo(stack_path);
% image_width_px = stack_info(1).Width;
% image_height_px = stack_info(1).Height;
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

% % Locate all the settings files in the dir
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

% Read in thunderstorm localisations
thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');

% Extract the required columns
thunderstorm_frame_array = thunderstorm_results.("frame");
thunderstorm_x_array = thunderstorm_results.("x [nm]");
thunderstorm_y_array = thunderstorm_results.("y [nm]");

% Loop over each image path and process

% % display ROI
% % Load the first frame from the stack
% first_frame = imread(stack_path, 1);
% first_frame = double(first_frame);
% 
% % Normalize the image for display
% display_image = (first_frame - min(first_frame(:))) / (max(first_frame(:)) - min(first_frame(:)));
% 
% % Convert ROI coordinates from nm to pixels directly (without using nm_to_px function)
% % Using manual conversion from center-based nm to top-left-based pixels
% image_center_x_px = image_width_px / 2;
% image_center_y_px = image_height_px / 2;
% 
% % Convert from nm (center origin) to pixels (top-left origin)
% ROI_min_x_px = image_center_x_px + (ROI_min_x / pixel_size_nm);
% ROI_max_x_px = image_center_x_px + (ROI_max_x / pixel_size_nm);
% ROI_min_y_px = image_center_y_px - (ROI_max_y / pixel_size_nm); % Note the negative sign for y
% ROI_max_y_px = image_center_y_px - (ROI_min_y / pixel_size_nm); % Note the negative sign for y
% 
% % Calculate rectangle position [left, top, width, height]
% ROI_left_px = min(ROI_min_x_px, ROI_max_x_px);
% ROI_top_px = min(ROI_min_y_px, ROI_max_y_px);
% ROI_width_px = abs(ROI_max_x_px - ROI_min_x_px);
% ROI_height_px = abs(ROI_max_y_px - ROI_min_y_px);
% 
% % Print out the values for debugging
% fprintf('ROI rectangle position: [%.2f, %.2f, %.2f, %.2f]\n', ...
%         ROI_left_px, ROI_top_px, ROI_width_px, ROI_height_px);
% fprintf('Image dimensions: %d x %d\n', image_width_px, image_height_px);
% 
% % Display the image
% figure;
% imshow(display_image);
% hold on;
% 
% % Draw the ROI rectangle in red with thicker line
% rectangle('Position', [ROI_left_px, ROI_top_px, ROI_width_px, ROI_height_px], ...
%           'EdgeColor', 'r', 'LineWidth', 3);
% 
% title('First Frame with ROI Highlighted');
% hold off;
% 
% pause(2);



% Loop over each frame in stack
for frame_index = 1:number_of_frames%length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
    % settings_path = settings_paths{frame_index};
    % 
    % % Read in parameters/settings if using simulated to get ground truth
    % evalc('run(settings_path)');
    % 
    % % filter only those within ROI
    % in_roi_mask = positionX_nm_array >= ROI_min_x & ...
    %           positionX_nm_array <= ROI_max_x & ...
    %           positionY_nm_array >= ROI_min_y & ...
    %           positionY_nm_array <= ROI_max_y;
    % positionX_nm_array = positionX_nm_array(in_roi_mask);
    % positionY_nm_array = positionY_nm_array(in_roi_mask);
    % angleInclination_array = angleInclination_array(in_roi_mask);
    % angleAzimuth_array = angleAzimuth_array(in_roi_mask);

    % Put hard-coded params in the right form
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveNA = objectiveNA;
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = patch_width_px;
    par.wavelength = Length(wavelength,'nm');
    par.refractiveIndices = refractiveIndices;
    par.pixelSensitivityMask = pixelSensitivityMask;
    par.nDiscretizationBFP = nDiscretizationBFP;
    par.backgroundNoise = backgroundNoise;
    par.nPhotons = nPhotons;

    % Load frame
    psf_image = imread(stack_path, frame_index);
    psf_image = double(psf_image); % Convert to double for calculations
    [image_height_px, image_width_px] = size(psf_image);
    
    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % Take thunderstorm results for that frame
    % Note: thunderstorm origin is top-left corner, so need to adjust to centre
    current_frame_mask = thunderstorm_frame_array == frame_index;
    current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
    current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

    % filter only those within ROI
    % convert thunderstorm coordinates to same reference frame as ROI
    % current_frame_x_array_centered = current_frame_x_array - image_width_nm/2;
    % current_frame_y_array_centered = -(current_frame_y_array - image_height_nm/2);
    in_roi_mask_ts = current_frame_x_array >= ROI_min_x & ...
                     current_frame_x_array <= ROI_max_x & ...
                     current_frame_y_array >= ROI_min_y & ...
                     current_frame_y_array <= ROI_max_y;
    current_frame_x_array = current_frame_x_array(in_roi_mask_ts);
    current_frame_y_array = current_frame_y_array(in_roi_mask_ts);

    % psf_image = psf_image ./ norm(psf_image);

    % outlined_image = psf_image;

    % Loop over each blob in a frame,
    % cropping out each blob
    % and running the fit on the cropped image

    for blob_index = 1:length(current_frame_x_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(current_frame_x_array));

        % Take thunderstorm localisation to be the patch centre
        patch_centre_x_nm = current_frame_x_array(blob_index);
        patch_centre_y_nm = current_frame_y_array(blob_index);

        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_width_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_height_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        % patch_width_px = floor(patch_width_nm / pixel_size_nm);
        % if mod(patch_width_px, 2) == 0
        %     patch_width_px = patch_width_px + 1; % Ensure odd number
        % end
        patch_width_px = patch_width_nm/pixel_size_nm;
        patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

        % Half width for calculations (integer division)
        half_width = floor(patch_width_px / 2);

        % Get image dimensions
        [image_height, image_width] = size(psf_image);

        % Adjust center coordinates to ensure the patch stays within image and keeps constant width
        adjusted_centre_x_px = min(max(half_width + 1, round(patch_centre_x_px)), image_width - half_width);
        adjusted_centre_y_px = min(max(half_width + 1, round(patch_centre_y_px)), image_height - half_width);

        % Calculate patch boundaries using the adjusted center
        patch_start_x_px = adjusted_centre_x_px - half_width;
        patch_start_y_px = adjusted_centre_y_px - half_width;
        patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        patch_end_y_px = patch_start_y_px + patch_width_px - 1;

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

        % Calculate expected position relative to patch center - beause of
        % the boundary checking stuff, it might not quite be zero anymore

        expected_offset_x = current_frame_x_array(blob_index) - actual_patch_center_x_nm;
        expected_offset_y = current_frame_y_array(blob_index) - actual_patch_center_y_nm;

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
        psfInit = PSF(par);
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
        parEst.parameterBounds.newangle1 = [-1, 1];
        parEst.parameterBounds.newangle2 = [-1, 1];
        parEst.parameterBounds.newangle3 = [0, 1];

        % Give centroid as initial position estimate
        parEst.parameterStartValues.x = Length(expected_offset_x, 'nm'); % Because by definition it will be in/near centre of patch
        parEst.parameterStartValues.y = Length(expected_offset_y, 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
        parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
        parEst.parameterStartValues.newangle3 = rand();

        % Now optimising over photon count too
        photon_number = par.nPhotons;
        photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
        parEst.parameterStartValues.photons = photon_estimate;
        parEst.parameterBounds.photons = [1, 1e10];%[photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_estimate = fitResult.estimatesPositionDefocus.ML(7);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_TS = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % Finding errors in the simulated data is a bit more complicated
        % because of thigns like multiple detections etc.
        % so for each localisation, find nearest ground truth to the
        % current patch and pick the smallest one. good enough for now.

        % % Find index of closest blob to patch centre and use that as ground truth
        % distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
        % [min_distance, min_index] = min(distances_to_ground_truth);

        % Use this nearest blob as ground truth for this localisation
        positionX_nm_TS = current_frame_x_array(blob_index);
        positionY_nm_TS = current_frame_y_array(blob_index);
        % angleInclination_TS = angleInclination_array(min_index);
        % angleAzimuth_TS = angleAzimuth_array(min_index);

        positionX_nm_error = positionX_nm_TS - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_TS - positionY_nm_estimate;
        % angleInclination_error = angleInclination_TS - angleInclination_estimate;
        % angleAzimuth_error = angleAzimuth_TS - angleAzimuth_estimate;
        % photons_error = photon_number - photons_estimate;
        % % objective_function_error = objective_function_TS - objective_function_estimate;

        % Append results for each blob to an array for this frame
        positionX_nm_TSs_frame(blob_index) = positionX_nm_TS;
        positionY_nm_TSs_frame(blob_index) = positionY_nm_TS;
        % angleInclination_TSs_frame(blob_index) = angleInclination_TS;
        % angleAzimuth_TSs_frame(blob_index) = angleAzimuth_TS;
        % photons_TSs_frame(blob_index) = photon_number;

        positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_estimate;

        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        % angleInclination_errors_frame(blob_index) = angleInclination_error;
        % angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        % photons_errors_frame(blob_index) = photons_error;
        % % objective_function_TSs_frame(blob_index) = 999999;%objective_function_TS;
        % % objective_function_estimates_frame(blob_index) = 999999;%objective_function_estimate;
        % % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;

        frames_array_frame(blob_index) = frame_index;

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    % Append results for each frame to an array for this whole thing

    positionX_nm_TSs{frame_index} = positionX_nm_TSs_frame;
    positionY_nm_TSs{frame_index} = positionY_nm_TSs_frame;
    % angleInclination_TSs{frame_index} = angleInclination_TSs_frame;
    % angleAzimuth_TSs{frame_index} = angleAzimuth_TSs_frame;
    % photons_TSs{frame_index} = photons_TSs_frame;%photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimates_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimates_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    % angleInclination_errors{frame_index} = angleInclination_errors_frame;
    % angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    % photons_errors{frame_index} = photons_errors_frame;

    % objective_function_TSs{frame_index} = objective_function_TSs_frame;
    % objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;

    frames_array{frame_index} = frames_array_frame;

    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % % Flatten the array (convert cell array to a numeric array)
    positionX_nm_TSs_flat = [positionX_nm_TSs{:}];
    positionY_nm_TSs_flat = [positionY_nm_TSs{:}];
    % angleInclination_TSs_flat = [angleInclination_TSs{:}];
    % angleAzimuth_TSs_flat = [angleAzimuth_TSs{:}];
    % photons_TSs_flat = [photons_TSs{:}];

    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    photons_estimates_flat = [photons_estimates{:}];

    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    % angleInclination_errors_flat = [angleInclination_errors{:}];
    % angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    % photons_errors_flat = [photons_errors{:}];

    % objective_function_TSs_flat = [objective_function_TSs{:}];
    % objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

    frames_array_flat = [frames_array{:}];

    % Save data about this stack of images (all frames compressed into one long array)

    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'frame = [%s]\n', sprintf('%i, ', frames_array_flat(1:end)));
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_TSs_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_TSs_flat(1:end)));
    % fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_TSs_flat(1:end)));
    % fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_TSs_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    % fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    % fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_TSs_flat(1:end)));
    % fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    % fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_TSs_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    % fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);







    % % Visualisation
    % 
    % figure;
    % subplot(1,2,1);
    % % Normalize the image for display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image);
    % hold on;
    % 
    % % Get coordinates for ground truth and estimated positions
    % gt_x_px = zeros(length(positionX_nm_TSs_frame), 1);
    % gt_y_px = zeros(length(positionY_nm_TSs_frame), 1);
    % est_x_px = zeros(length(positionX_nm_estimates_frame), 1);
    % est_y_px = zeros(length(positionY_nm_estimates_frame), 1);
    % 
    % for i = 1:length(positionX_nm_TSs_frame)
    %     % Convert from nm to pixels
    %     gt_x_px(i) = nm_to_px(positionX_nm_TSs_frame(i), pixel_size_nm, image_width_px, 'x');
    %     gt_y_px(i) = nm_to_px(positionY_nm_TSs_frame(i), pixel_size_nm, image_height_px, 'y');
    %     est_x_px(i) = nm_to_px(positionX_nm_estimates_frame(i), pixel_size_nm, image_width_px, 'x');
    %     est_y_px(i) = nm_to_px(positionY_nm_estimates_frame(i), pixel_size_nm, image_height_px, 'y');
    % 
    %     % Draw line between ground truth and estimate
    %     plot([gt_x_px(i), est_x_px(i)], [gt_y_px(i), est_y_px(i)], 'y-', 'LineWidth', 1);
    % end
    % 
    % % Plot ground truth positions (green)
    % plot(gt_x_px, gt_y_px, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    % % Plot estimated positions (red)
    % plot(est_x_px, est_y_px, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    % 
    % title(sprintf('Frame %d: Ground Truth (o) vs Estimates (+)', frame_index));
    % legend('Error Lines', 'Ground Truth', 'Estimates');
    % 
    % % Subplot for error histogram
    % subplot(1,2,2);
    % % Calculate position errors in nm
    % pos_errors = sqrt(positionX_nm_errors_frame.^2 + positionY_nm_errors_frame.^2);
    % histogram(pos_errors, 20);
    % title('Position Error Histogram (nm)');
    % xlabel('Error (nm)');
    % ylabel('Frequency');
    % xlim([0, max(pos_errors) + 10]);
    % 
    % % Add some statistics to the plot
    % mean_err = mean(pos_errors);
    % median_err = median(pos_errors);
    % max_err = max(pos_errors);
    % text(0.05, 0.9, sprintf('Mean Error: %.2f nm', mean_err), 'Units', 'normalized');
    % text(0.05, 0.85, sprintf('Median Error: %.2f nm', median_err), 'Units', 'normalized');
    % text(0.05, 0.8, sprintf('Max Error: %.2f nm', max_err), 'Units', 'normalized');
    % 
    % % Make the figure bigger for better visibility
    % set(gcf, 'Position', [100, 100, 1200, 500]);
    % 
    % % Optional: save the figure if desired
    % saveas(gcf, sprintf('%s/visualization_frame_%d.png', stack_dir, frame_index));
    % 
    % % Pause briefly to allow viewing before moving to next frame
    % pause(0.5);







end % end loop over frames



