
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
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_vary_az_reduced/sims_highN/all/';

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_vary_az_reduced/fitting_results_hinterer_highN_swarm_reparam.py';




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

% Initialise overall results arrays
positionX_true = [];
positionY_true = [];
angleInclination_true = [];
angleAzimuth_true = [];

positionX_estimates = [];
positionY_estimates = [];
angleInclination_estimates = [];
angleAzimuth_estimates = [];

positionX_errors = [];
positionY_errors = [];
angleInclination_errors = [];
angleAzimuth_errors = [];

objective_function_trues = [];
objective_function_estimates = [];
objective_function_errors = [];

% Loop over each frame
for frame_index = 1:length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, length(frame_paths));
    frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in ground truth data:
    % run(settings_path);
    evalc('run(settings_path)');

    par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = image_size_px;
    par.wavelength = Length(wavelength,'nm');
    patch_width_nm = 1000; % size of patch around blob to consider
    patch_width_px = patch_width_nm/pixel_size_nm;

    % And convert to image coords for drawing later
    positionX_px_array = [];
    positionY_px_array = [];
    for i=1:number_of_spots
        positionX_px_array(end+1) = nm_to_px(positionX_nm_array(i), pixel_size_nm, image_size_px, 'x');
        positionY_px_array(end+1) = nm_to_px(positionY_nm_array(i), pixel_size_nm, image_size_px, 'y');
    end

    f_array = ones(1, number_of_spots);
    x_array = positionX_px_array; % replace this with output of some blob-finding algorithm
    y_array = positionY_px_array;

    % Load frame
    psf_image = imread(frame_path);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    psf_image = psf_image ./ norm(psf_image);

    outlined_image = psf_image;

    % Loop over each blob in a frame,
    % masking off everything but the blob
    % and running the fit on the masked image

    % Initialise results arrays
    positionX_nm_estimate_frame = [];
    positionY_nm_estimate_frame = [];
    angleInclination_estimates_frame = [];
    angleAzimuth_estimates_frame = [];
    positionX_nm_errors_frame = [];
    positionY_nm_errors_frame = [];
    angleInclination_errors_frame = [];
    angleAzimuth_errors_frame = [];

    for blob_index = 1:length(x_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(x_array));

        centroid_x_nm = positionX_nm_array(blob_index) + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        centroid_y_nm = positionY_nm_array(blob_index) + (rand - 0.5) * 200;

        % % If masking, use this
        % centroid_x_px = positionX_px_array(blob_index);
        % centroid_y_px = positionY_px_array(blob_index);
        % 
        % patch_start_x_nm = centroid_x_nm - patch_width_nm/2;
        % patch_start_y_nm = centroid_y_nm - patch_width_nm/2;
        % patch_end_x_nm = centroid_x_nm + patch_width_nm/2;
        % patch_end_y_nm = centroid_y_nm + patch_width_nm/2;
        % 
        % patch_start_x_px = round(centroid_x_px - (patch_width_px + 1)/2);
        % patch_start_y_px = round(centroid_y_px - (patch_width_px + 1)/2);
        % patch_end_x_px = round(centroid_x_px + (patch_width_px + 1)/2);
        % patch_end_y_px = round(centroid_y_px + (patch_width_px + 1)/2);
        % 
        % % Generate throwaway PSF object, later replace the image with our masked image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % 
        % % Create masked image, where everything outside patch = 0
        % masked_psf_image = psf_image;
        % [image_height, image_width] = size(masked_psf_image);  % Get the image size
        % 
        % mask_x = (1:image_width) >= patch_start_x_px & (1:image_width) <= patch_end_x_px;  % x indices within the patch
        % mask_y = (1:image_height) >= patch_start_y_px & (1:image_height) <= patch_end_y_px;  % y indices within the patch
        % 
        % masked_psf_image(~mask_y, :) = 0;
        % masked_psf_image(:, ~mask_x) = 0;
        % 
        % % Replace image in PSF object with masked image
        % psfInit.image = masked_psf_image;
       
        % If not masking, use this
        % Generate throwaway PSF object, replace the image with our image
        par.position = Length([0 0 0], 'nm');
        psfInit = PSF(par);
        psfInit.image = psf_image;

        % % Give patch as localisation bounds
        % parEst.parameterBounds.x = Length([patch_start_x_nm patch_end_x_nm], 'nm');
        % parEst.parameterBounds.y = Length([patch_start_y_nm patch_end_y_nm], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % parEst.parameterBounds.azimuth = [-0.0174533 (2*pi)+0.0174533];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2; % add angle optimiser
        % parEst.parameterStartValues.azimuth = rand() * 2*pi; % add angle optimiser

        % % Run the fit - no angle optimising
        % parEst.angleInclinationEstimate = 0;%angleInclination_array(blob_index); % given orientations in this case. set to 0 otherwise
        % parEst.angleAzimuthEstimate = 0;%angleAzimuth_array(blob_index); % given orientations in this case. set to 0 otherwise
        % tic;
        % fitResult = FitPSF(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate;

        % % Run the fit - Gaussian model, no angle optimising
        % parEst.angleInclinationEstimate = pi/4;%angleInclination_array(blob_index); % given orientations in this case. set to 0 otherwise
        % parEst.angleAzimuthEstimate = pi/4;%angleAzimuth_array(blob_index); % given orientations in this case. set to 0 otherwise
        % tic;
        % fitResult = FitPSF_gaussian(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate;

        % Run the fit - with angle optimising
        % tic;
        % % fitResult = FitPSF_angles(psfInit, parEst);
        % % fitResult = FitPSF_angles_global(psfInit, parEst);
        % % fitResult = FitPSF_gaussian(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);

        % % Run the fit - with angle optimising, not polar, fminsearch
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % tic;
        % fitResult = FitPSF_notpolar_global(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(4);
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);




        % % Run the fit - given nothing, fmincon
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % parEst.parameterBounds.azimuth = [0 2*pi];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2;
        % parEst.parameterStartValues.azimuth = rand() * 2*pi;
        % % tic;
        % fitResult = FitPSF_ML(psfInit, parEst);
        % % fitResult = FitPSF_PS_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.ML(6);

        % % Run the fit - given nothing, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-1 1], 'nm');
        % parEst.parameterBounds.inclination = [0, pi/2];
        % parEst.parameterBounds.azimuth = [0, 2*pi];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2;
        % parEst.parameterStartValues.azimuth = rand() * 2 * pi;
        % % tic;
        % fitResult = FitPSF_PS(psfInit, parEst);
        % % fitResult = FitPSF_PS_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = mod(fitResult.estimatesPositionDefocus.PS(4), pi/2);
        % angleAzimuth_estimate = mod(fitResult.estimatesPositionDefocus.PS(5), 2*pi);
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(6);



        % Run the fit - given nothing, particle swarm - reparameterised
        parEst.parameterBounds.x = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.sinInclination = [-1, 1];
        parEst.parameterBounds.cosAzimuth = [-1, 1];
        parEst.parameterBounds.sinAzimuth = [-1, 1];
        % Give centroid as initial position estimate
        parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.sinInclination = 2*rand()-0.5;
        parEst.parameterStartValues.cosAzimuth = 2*(rand()-0.5);
        parEst.parameterStartValues.sinAzimuth = 2*(rand()-0.5);
        % tic;
        fitResult = FitPSF_PS_reparam(psfInit, parEst);
        % fitResult = FitPSF_PS_gaussian(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = 0.5*asin(fitResult.estimatesPositionDefocus.PS(4));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.PS(6), fitResult.estimatesPositionDefocus.PS(5));
        positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        objective_function_estimate = fitResult.estimatesPositionDefocus.PS(7);



        % % Run the fit - given polar, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.azimuth = [0 2*pi];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.azimuth = rand() * 2*pi;
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % % tic;
        % fitResult = FitPSF_PS_givenpolar(psfInit, parEst);
        % % fitResult = FitPSF_PS_givenpolar_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.PS(4); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(5);

        % % Run the fit - given azimuth, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2;
        % parEst.angleAzimuthEstimate = angleAzimuth_array(blob_index); % given orientations in this case
        % % tic;
        % fitResult = FitPSF_PS_givenazimuth(psfInit, parEst);
        % % fitResult = FitPSF_PS_givenpolar_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate;
        % angleInclination_estimate = fitResult.estimatesPositionDefocus.PS(4); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(5);

        % % Run the fit - given orientation, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % parEst.angleAzimuthEstimate = angleAzimuth_array(blob_index); % given orientations in this case
        % % tic;
        % fitResult = FitPSF_PS_givenorientation(psfInit, parEst);
        % % fitResult = FitPSF_PS_givenorientation_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate; 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(4);

        % Now we need the value of the objective function at ground truth
        ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % % Run the fit - with angle optimising, not polar, particle swarm,
        % % optimise sin(az) cos(az), not just az
        % parEst.parameterBounds.x = Length([patch_start_x_nm patch_end_x_nm], 'nm');
        % parEst.parameterBounds.y = Length([patch_start_y_nm patch_end_y_nm], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % parEst.parameterBounds.sin_azimuth = [-1 1];
        % parEst.parameterBounds.cos_azimuth = [-1 1];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = 0; % add angle optimiser
        % parEst.parameterStartValues.sin_azimuth = 0; % add angle optimiser
        % parEst.parameterStartValues.cos_azimuth = 0; % add angle optimiser
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % tic;
        % fitResult = FitPSF_PS_sine(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.PS(4); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);

        % wrap these in abs if doing that again
        positionX_nm_error = positionX_nm_array(blob_index) - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_array(blob_index) - positionY_nm_estimate;
        angleInclination_error = angleInclination_array(blob_index) - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_array(blob_index) - angleAzimuth_estimate;
        objective_function_error = objective_function_true - objective_function_estimate;

        % Store results
        positionX_nm_estimate_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimate_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        objective_function_trues_frame(blob_index) = objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        objective_function_errors_frame(blob_index) = objective_function_error;

        % fprintf('    Error: (Δx, Δy) = (%i, %i)\n', round(positionX_nm_error), round(positionY_nm_error));
        % fprintf('    Error: (Δθ, Δφ) = (%i°, %i°)\n', round(angleInclination_error*180/pi), round(angleAzimuth_error*180/pi));
        % 
        % % Draw a red rectangle on the original image for the patch
        % current_color = 'red';
        % outlined_image = insertShape(outlined_image, 'Rectangle', ...
        %                              [patch_start_x_px, patch_start_y_px, patch_width_px, patch_width_px], ...
        %                              'Color', current_color, 'LineWidth', 1);
        % outlined_image = insertShape(outlined_image, 'Circle', [centroid_x_px, centroid_y_px, 0], 'Color', 'red', 'LineWidth', 1);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    positionX_nm_true{frame_index} = positionX_nm_array;
    positionY_nm_true{frame_index} = positionY_nm_array;
    angleInclination_true{frame_index} = angleInclination_array;
    angleAzimuth_true{frame_index} = angleAzimuth_array;

    positionX_nm_estimates{frame_index} = positionX_nm_estimate_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimate_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;

    objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    objective_function_errors{frame_index} = objective_function_errors_frame;
        
    
    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_true_flat = [positionX_nm_true{:}];
    positionY_nm_true_flat = [positionY_nm_true{:}];
    angleInclination_true_flat = [angleInclination_true{:}];
    angleAzimuth_true_flat = [angleAzimuth_true{:}];
    
    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    
    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    
    objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    objective_function_errors_flat = [objective_function_errors{:}];
    
    % Save data about this stack of images (all frames compressed into one long array)
    
    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_true_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_true_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_true_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_true_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames









%% Low photon number


%% Fitting multiple PSFs in a single frame

% Same as loop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;


%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_vary_az_reduced/sims/all/';

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_vary_az_reduced/fitting_results_hinterer_lowN_swarm_reparam.py';




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

% Initialise overall results arrays
positionX_true = [];
positionY_true = [];
angleInclination_true = [];
angleAzimuth_true = [];

positionX_estimates = [];
positionY_estimates = [];
angleInclination_estimates = [];
angleAzimuth_estimates = [];

positionX_errors = [];
positionY_errors = [];
angleInclination_errors = [];
angleAzimuth_errors = [];

objective_function_trues = [];
objective_function_estimates = [];
objective_function_errors = [];

% Loop over each frame
for frame_index = 1:length(frame_paths)%randperm(length(frame_paths), 40)

    tic; % timing each frame

    fprintf('----------\n');
    fprintf('FRAME %d/%d\n', frame_index, length(frame_paths));
    frame_path = frame_paths{frame_index};
    settings_path = settings_paths{frame_index};

    % Read in ground truth data:
    % run(settings_path);
    evalc('run(settings_path)');

    par.pixelSensitivityMask = PixelSensitivity.uniform(9);
    par.pixelSize = Length(pixel_size_nm,'nm');
    par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
    par.nPixels = image_size_px;
    par.wavelength = Length(wavelength,'nm');
    patch_width_nm = 1000; % size of patch around blob to consider
    patch_width_px = patch_width_nm/pixel_size_nm;

    % And convert to image coords for drawing later
    positionX_px_array = [];
    positionY_px_array = [];
    for i=1:number_of_spots
        positionX_px_array(end+1) = nm_to_px(positionX_nm_array(i), pixel_size_nm, image_size_px, 'x');
        positionY_px_array(end+1) = nm_to_px(positionY_nm_array(i), pixel_size_nm, image_size_px, 'y');
    end

    f_array = ones(1, number_of_spots);
    x_array = positionX_px_array; % replace this with output of some blob-finding algorithm
    y_array = positionY_px_array;

    % Load frame
    psf_image = imread(frame_path);
    psf_image = double(psf_image); % Convert to double for calculations

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    psf_image = psf_image ./ norm(psf_image);

    outlined_image = psf_image;

    % Loop over each blob in a frame,
    % masking off everything but the blob
    % and running the fit on the masked image

    % Initialise results arrays
    positionX_nm_estimate_frame = [];
    positionY_nm_estimate_frame = [];
    angleInclination_estimates_frame = [];
    angleAzimuth_estimates_frame = [];
    positionX_nm_errors_frame = [];
    positionY_nm_errors_frame = [];
    angleInclination_errors_frame = [];
    angleAzimuth_errors_frame = [];

    for blob_index = 1:length(x_array)

        fprintf(' ∘ Blob %d/%d\n', blob_index, length(x_array));

        centroid_x_nm = positionX_nm_array(blob_index) + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        centroid_y_nm = positionY_nm_array(blob_index) + (rand - 0.5) * 200;

        % % If masking, use this
        % centroid_x_px = positionX_px_array(blob_index);
        % centroid_y_px = positionY_px_array(blob_index);
        % 
        % patch_start_x_nm = centroid_x_nm - patch_width_nm/2;
        % patch_start_y_nm = centroid_y_nm - patch_width_nm/2;
        % patch_end_x_nm = centroid_x_nm + patch_width_nm/2;
        % patch_end_y_nm = centroid_y_nm + patch_width_nm/2;
        % 
        % patch_start_x_px = round(centroid_x_px - (patch_width_px + 1)/2);
        % patch_start_y_px = round(centroid_y_px - (patch_width_px + 1)/2);
        % patch_end_x_px = round(centroid_x_px + (patch_width_px + 1)/2);
        % patch_end_y_px = round(centroid_y_px + (patch_width_px + 1)/2);
        % 
        % % Generate throwaway PSF object, later replace the image with our masked image
        % par.position = Length([0 0 0], 'nm');
        % psfInit = PSF(par);
        % 
        % % Create masked image, where everything outside patch = 0
        % masked_psf_image = psf_image;
        % [image_height, image_width] = size(masked_psf_image);  % Get the image size
        % 
        % mask_x = (1:image_width) >= patch_start_x_px & (1:image_width) <= patch_end_x_px;  % x indices within the patch
        % mask_y = (1:image_height) >= patch_start_y_px & (1:image_height) <= patch_end_y_px;  % y indices within the patch
        % 
        % masked_psf_image(~mask_y, :) = 0;
        % masked_psf_image(:, ~mask_x) = 0;
        % 
        % % Replace image in PSF object with masked image
        % psfInit.image = masked_psf_image;
       
        % If not masking, use this
        % Generate throwaway PSF object, replace the image with our image
        par.position = Length([0 0 0], 'nm');
        psfInit = PSF(par);
        psfInit.image = psf_image;

        % % Give patch as localisation bounds
        % parEst.parameterBounds.x = Length([patch_start_x_nm patch_end_x_nm], 'nm');
        % parEst.parameterBounds.y = Length([patch_start_y_nm patch_end_y_nm], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % parEst.parameterBounds.azimuth = [-0.0174533 (2*pi)+0.0174533];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2; % add angle optimiser
        % parEst.parameterStartValues.azimuth = rand() * 2*pi; % add angle optimiser

        % % Run the fit - no angle optimising
        % parEst.angleInclinationEstimate = 0;%angleInclination_array(blob_index); % given orientations in this case. set to 0 otherwise
        % parEst.angleAzimuthEstimate = 0;%angleAzimuth_array(blob_index); % given orientations in this case. set to 0 otherwise
        % tic;
        % fitResult = FitPSF(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate;

        % % Run the fit - Gaussian model, no angle optimising
        % parEst.angleInclinationEstimate = pi/4;%angleInclination_array(blob_index); % given orientations in this case. set to 0 otherwise
        % parEst.angleAzimuthEstimate = pi/4;%angleAzimuth_array(blob_index); % given orientations in this case. set to 0 otherwise
        % tic;
        % fitResult = FitPSF_gaussian(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate;

        % Run the fit - with angle optimising
        % tic;
        % % fitResult = FitPSF_angles(psfInit, parEst);
        % % fitResult = FitPSF_angles_global(psfInit, parEst);
        % % fitResult = FitPSF_gaussian(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);

        % % Run the fit - with angle optimising, not polar, fminsearch
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % tic;
        % fitResult = FitPSF_notpolar_global(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(4);
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);




        % % Run the fit - given nothing, fmincon
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % parEst.parameterBounds.azimuth = [0 2*pi];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2;
        % parEst.parameterStartValues.azimuth = rand() * 2*pi;
        % % tic;
        % fitResult = FitPSF_ML(psfInit, parEst);
        % % fitResult = FitPSF_PS_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.ML(6);

        % % Run the fit - given nothing, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-1 1], 'nm');
        % parEst.parameterBounds.inclination = [0, pi/2];
        % parEst.parameterBounds.azimuth = [0, 2*pi];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2;
        % parEst.parameterStartValues.azimuth = rand() * 2 * pi;
        % % tic;
        % fitResult = FitPSF_PS(psfInit, parEst);
        % % fitResult = FitPSF_PS_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = mod(fitResult.estimatesPositionDefocus.PS(4), pi/2);
        % angleAzimuth_estimate = mod(fitResult.estimatesPositionDefocus.PS(5), 2*pi);
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(6);



        % Run the fit - given nothing, particle swarm - reparameterised
        parEst.parameterBounds.x = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.sinInclination = [-1, 1];
        parEst.parameterBounds.cosAzimuth = [-1, 1];
        parEst.parameterBounds.sinAzimuth = [-1, 1];
        % Give centroid as initial position estimate
        parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.sinInclination = 2*rand()-0.5;
        parEst.parameterStartValues.cosAzimuth = 2*(rand()-0.5);
        parEst.parameterStartValues.sinAzimuth = 2*(rand()-0.5);
        % tic;
        fitResult = FitPSF_PS_reparam(psfInit, parEst);
        % fitResult = FitPSF_PS_gaussian(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = 0.5*asin(fitResult.estimatesPositionDefocus.PS(4));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.PS(6), fitResult.estimatesPositionDefocus.PS(5));
        positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        objective_function_estimate = fitResult.estimatesPositionDefocus.PS(7);



        % % Run the fit - given polar, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.azimuth = [0 2*pi];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.azimuth = rand() * 2*pi;
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % % tic;
        % fitResult = FitPSF_PS_givenpolar(psfInit, parEst);
        % % fitResult = FitPSF_PS_givenpolar_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.PS(4); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(5);

        % % Run the fit - given azimuth, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = rand() * pi/2;
        % parEst.angleAzimuthEstimate = angleAzimuth_array(blob_index); % given orientations in this case
        % % tic;
        % fitResult = FitPSF_PS_givenazimuth(psfInit, parEst);
        % % fitResult = FitPSF_PS_givenpolar_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate;
        % angleInclination_estimate = fitResult.estimatesPositionDefocus.PS(4); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(5);

        % % Run the fit - given orientation, particle swarm
        % parEst.parameterBounds.x = Length([centroid_x_nm-200 centroid_x_nm+200], 'nm');
        % parEst.parameterBounds.y = Length([centroid_y_nm-200 centroid_y_nm+200], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % parEst.angleAzimuthEstimate = angleAzimuth_array(blob_index); % given orientations in this case
        % % tic;
        % fitResult = FitPSF_PS_givenorientation(psfInit, parEst);
        % % fitResult = FitPSF_PS_givenorientation_gaussian(psfInit, parEst);
        % % elapsed_time = toc;
        % % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = parEst.angleAzimuthEstimate; 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);
        % objective_function_estimate = fitResult.estimatesPositionDefocus.PS(4);

        % Now we need the value of the objective function at ground truth
        ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % % Run the fit - with angle optimising, not polar, particle swarm,
        % % optimise sin(az) cos(az), not just az
        % parEst.parameterBounds.x = Length([patch_start_x_nm patch_end_x_nm], 'nm');
        % parEst.parameterBounds.y = Length([patch_start_y_nm patch_end_y_nm], 'nm');
        % parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        % parEst.parameterBounds.inclination = [0 pi/2];
        % parEst.parameterBounds.sin_azimuth = [-1 1];
        % parEst.parameterBounds.cos_azimuth = [-1 1];
        % % Give centroid as initial position estimate
        % parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        % parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        % parEst.parameterStartValues.defocus = Length(0, 'nm');
        % parEst.parameterStartValues.inclination = 0; % add angle optimiser
        % parEst.parameterStartValues.sin_azimuth = 0; % add angle optimiser
        % parEst.parameterStartValues.cos_azimuth = 0; % add angle optimiser
        % parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
        % tic;
        % fitResult = FitPSF_PS_sine(psfInit, parEst);
        % elapsed_time = toc;
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        % angleInclination_estimate = parEst.angleInclinationEstimate;
        % angleAzimuth_estimate = fitResult.estimatesPositionDefocus.PS(4); 
        % positionX_nm_estimate = fitResult.estimatesPositionDefocus.PS(1);
        % positionY_nm_estimate = fitResult.estimatesPositionDefocus.PS(2);

        % wrap these in abs if doing that again
        positionX_nm_error = positionX_nm_array(blob_index) - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_array(blob_index) - positionY_nm_estimate;
        angleInclination_error = angleInclination_array(blob_index) - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_array(blob_index) - angleAzimuth_estimate;
        objective_function_error = objective_function_true - objective_function_estimate;

        % Store results
        positionX_nm_estimate_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimate_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        objective_function_trues_frame(blob_index) = objective_function_true;
        objective_function_estimates_frame(blob_index) = objective_function_estimate;
        objective_function_errors_frame(blob_index) = objective_function_error;

        % fprintf('    Error: (Δx, Δy) = (%i, %i)\n', round(positionX_nm_error), round(positionY_nm_error));
        % fprintf('    Error: (Δθ, Δφ) = (%i°, %i°)\n', round(angleInclination_error*180/pi), round(angleAzimuth_error*180/pi));
        % 
        % % Draw a red rectangle on the original image for the patch
        % current_color = 'red';
        % outlined_image = insertShape(outlined_image, 'Rectangle', ...
        %                              [patch_start_x_px, patch_start_y_px, patch_width_px, patch_width_px], ...
        %                              'Color', current_color, 'LineWidth', 1);
        % outlined_image = insertShape(outlined_image, 'Circle', [centroid_x_px, centroid_y_px, 0], 'Color', 'red', 'LineWidth', 1);

    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    positionX_nm_true{frame_index} = positionX_nm_array;
    positionY_nm_true{frame_index} = positionY_nm_array;
    angleInclination_true{frame_index} = angleInclination_array;
    angleAzimuth_true{frame_index} = angleAzimuth_array;

    positionX_nm_estimates{frame_index} = positionX_nm_estimate_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimate_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;

    objective_function_trues{frame_index} = objective_function_trues_frame;
    objective_function_estimates{frame_index} = objective_function_estimates_frame;
    objective_function_errors{frame_index} = objective_function_errors_frame;
        
    
    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_true_flat = [positionX_nm_true{:}];
    positionY_nm_true_flat = [positionY_nm_true{:}];
    angleInclination_true_flat = [angleInclination_true{:}];
    angleAzimuth_true_flat = [angleAzimuth_true{:}];
    
    positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
    positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
    angleInclination_estimates_flat = [angleInclination_estimates{:}];
    angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
    
    positionX_nm_errors_flat = [positionX_nm_errors{:}];
    positionY_nm_errors_flat = [positionY_nm_errors{:}];
    angleInclination_errors_flat = [angleInclination_errors{:}];
    angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
    
    objective_function_trues_flat = [objective_function_trues{:}];
    objective_function_estimates_flat = [objective_function_estimates{:}];
    objective_function_errors_flat = [objective_function_errors{:}];
    
    % Save data about this stack of images (all frames compressed into one long array)
    
    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_true_flat(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_true_flat(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_true_flat(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_true_flat(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
    fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames


