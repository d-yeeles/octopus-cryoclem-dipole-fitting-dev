% 
% %% Fitting multiple PSFs in a single frame
% 
% % Same as loop-test/fit_multiple.m, but feeding it the known simulated
% % positions rather than using thunderstorm to estimate them.
% 
% close all;
% clear all;
% 
% addpath(genpath('../'));
% 
% %% ----------
% %% Fit
% %% ----------
% 
% % Input params
% frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/sims_1spot/massive/';
% 
% % Output file
% fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/results_1spot_old_massive.py';
% 
% % Model to use
% model = 'hinterer';
% 
% 
% 
% % Function to create the current PSF (same as `createFitPSF` method in the class)
% % function currentPSF = createFitPSF(psfEstimate, params)
% %     psfEstimate.position = Length([params(1:2), 0], 'nm');
% %     psfEstimate.defocus = Length(params(3), 'nm');
% %     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
% %     noiseEstimate = psfEstimate.backgroundNoise;
% %     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% % 
% %     % Simulate the PSF
% %     bfp = BackFocalPlane(psfEstimate);
% %     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
% %     psfEstimate.backFocalPlane = bfp;
% % 
% %     % Apply phase mask
% %     psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
% %     psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
% % 
% %     % Apply attenuation mask
% %     psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
% %     psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
% % 
% %     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% % 
% %     for k=1:size(psfEstimate.stageDrift.motion,1)
% %         % Apply aberrations
% %         aberrations = getAberrations(psfEstimate,k);
% %         aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
% % 
% %         % Get image from BFP field
% %         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
% %     end
% % 
% %     currentPsf = adjustExcitation(psfEstimate, currentPsf);
% %     currentPsf = applyShotNoise(psfEstimate, currentPsf);
% %     currentPsf = addBackgroundNoise(psfEstimate, currentPsf);
% % 
% %     totalIntensity = sum(currentPsf,'all');
% %     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
% %     currentFitPSF = currentPsf ./ norm(currentPsf);
% % 
% %     currentPSF = currentFitPSF;
% % end
% % function currentPSF = createFitPSF(psfEstimate, params)
% % 
% %     psfEstimate.position = Length([params(1:2), 0], 'nm');
% %     psfEstimate.defocus = Length(params(3), 'nm');
% %     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
% %     noiseEstimate = psfEstimate.backgroundNoise;
% %     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% % 
% %     bfp = BackFocalPlane(psfEstimate);
% %     % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
% %     psfEstimate.backFocalPlane = bfp;
% % 
% %     % dave
% %     psfEstimate.fieldBFP.x = bfp.electricField.x;
% %     psfEstimate.fieldBFP.y = bfp.electricField.y;
% % 
% %     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% %     for k=1:size(psfEstimate.stageDrift.motion,1)
% %         aberrationCoeffs = getAberrations(psfEstimate,k);
% %         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
% %         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
% %     end
% % 
% %     totalIntensity = sum(currentPsf,'all');
% %     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
% %     currentFitPSF = currentPsf ./ norm(currentPsf);
% % 
% %     currentPSF = currentFitPSF;
% % 
% % end
% % 
% % function currentPSF = createFitPSF_gaussian(psfEstimate, params)
% % 
% %     psfEstimate.position = Length([params(1:2), 0], 'nm');
% %     psfEstimate.defocus = Length(params(3), 'nm');
% %     psfEstimate.dipole = Dipole(params(4), params(5));  % Use inclination
% %     noiseEstimate = psfEstimate.backgroundNoise;
% %     nPhotonEstimate = psfEstimate.nPhotons;%round(sum(sum(psfEstimate.image - noiseEstimate)));
% % 
% %     % bfp = BackFocalPlane(psfEstimate);
% %     bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
% %     psfEstimate.backFocalPlane = bfp;
% % 
% %     % dave
% %     psfEstimate.fieldBFP.x = bfp.electricField.x;
% %     psfEstimate.fieldBFP.y = bfp.electricField.y;
% % 
% %     currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
% %     for k=1:size(psfEstimate.stageDrift.motion,1)
% %         aberrationCoeffs = getAberrations(psfEstimate,k);
% %         fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
% %         currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
% %     end
% % 
% %     totalIntensity = sum(currentPsf,'all');
% %     currentPsf = currentPsf ./ totalIntensity * nPhotonEstimate + noiseEstimate;
% %     currentFitPSF = currentPsf ./ norm(currentPsf);
% % 
% %     currentPSF = currentFitPSF;
% % 
% % end
% 
% 
% 
% 
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
% 
% % Create the centroids in image coordinates
% % 'Read in' ground truth positions (this will be rpelaced by like reading in thunderstorm estimate or something)
% 
% % Loop over each image path and process
% 
% % Initialise overall results arrays
% positionX_true = [];
% positionY_true = [];
% angleInclination_true = [];
% angleAzimuth_true = [];
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
% 
% % Loop over each frame
% for frame_index = 1:length(frame_paths)%randperm(length(frame_paths), 40)
% 
%     tic; % timing each frame
% 
%     fprintf('----------\n');
%     fprintf('FRAME %d/%d\n', frame_index, length(frame_paths));
%     frame_path = frame_paths{frame_index};
%     settings_path = settings_paths{frame_index};
% 
%     % Read in ground truth data:
%     % run(settings_path);
%     evalc('run(settings_path)');
% 
%     par.pixelSensitivityMask = PixelSensitivity.uniform(9);
%     par.pixelSize = Length(pixel_size_nm,'nm');
%     par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
%     par.nPixels = image_size_px;
%     par.wavelength = Length(wavelength,'nm');
%     patch_width_nm = image_size_nm; % size of patch around blob to consider
%     patch_width_px = patch_width_nm/pixel_size_nm;
% 
%     % And convert to image coords for drawing later
%     positionX_px_array = [];
%     positionY_px_array = [];
%     for i=1:number_of_spots
%         positionX_px_array(end+1) = nm_to_px(positionX_nm_array(i), pixel_size_nm, image_size_px, 'x');
%         positionY_px_array(end+1) = nm_to_px(positionY_nm_array(i), pixel_size_nm, image_size_px, 'y');
%     end
% 
%     f_array = ones(1, number_of_spots);
%     x_array = positionX_px_array; % replace this with output of some blob-finding algorithm
%     y_array = positionY_px_array;
% 
%     % Load frame
%     psf_image = imread(frame_path); % if tif
%     psf_image = double(psf_image); % Convert to double for calculations if tif
%     % psf_image = readmatrix(frame_path); % if csv
% 
%     % % Clip values just if want to display
%     % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
%     % imshow(display_image)
% 
%     % psf_image = psf_image ./ norm(psf_image);
% 
%     outlined_image = psf_image;
% 
%     % Loop over each blob in a frame,
%     % masking off everything but the blob
%     % and running the fit on the masked image
% 
%     % Initialise results arrays
%     positionX_nm_estimate_frame = [];
%     positionY_nm_estimate_frame = [];
%     angleInclination_estimates_frame = [];
%     angleAzimuth_estimates_frame = [];
%     positionX_nm_errors_frame = [];
%     positionY_nm_errors_frame = [];
%     angleInclination_errors_frame = [];
%     angleAzimuth_errors_frame = [];
% 
%     for blob_index = 1:length(x_array)
% 
%         fprintf(' ∘ Blob %d/%d\n', blob_index, length(x_array));
% 
%         centroid_x_nm = positionX_nm_array(blob_index) + (rand - 0.5) * 200; % this is simulating some other blob-finding result
%         centroid_y_nm = positionY_nm_array(blob_index) + (rand - 0.5) * 200;
% 
%         % % If masking, use this
%         % centroid_x_px = positionX_px_array(blob_index);
%         % centroid_y_px = positionY_px_array(blob_index);
%         % 
%         % patch_start_x_nm = centroid_x_nm - patch_width_nm/2;
%         % patch_start_y_nm = centroid_y_nm - patch_width_nm/2;
%         % patch_end_x_nm = centroid_x_nm + patch_width_nm/2;
%         % patch_end_y_nm = centroid_y_nm + patch_width_nm/2;
%         % 
%         % patch_start_x_px = round(centroid_x_px - (patch_width_px + 1)/2);
%         % patch_start_y_px = round(centroid_y_px - (patch_width_px + 1)/2);
%         % patch_end_x_px = round(centroid_x_px + (patch_width_px + 1)/2);
%         % patch_end_y_px = round(centroid_y_px + (patch_width_px + 1)/2);
%         % 
%         % % Generate throwaway PSF object, later replace the image with our masked image
%         % par.position = Length([0 0 0], 'nm');
%         % psfInit = PSF(par);
%         % 
%         % % Create masked image, where everything outside patch = 0
%         % masked_psf_image = psf_image;
%         % [image_height, image_width] = size(masked_psf_image);  % Get the image size
%         % 
%         % mask_x = (1:image_width) >= patch_start_x_px & (1:image_width) <= patch_end_x_px;  % x indices within the patch
%         % mask_y = (1:image_height) >= patch_start_y_px & (1:image_height) <= patch_end_y_px;  % y indices within the patch
%         % 
%         % masked_psf_image(~mask_y, :) = 0;
%         % masked_psf_image(:, ~mask_x) = 0;
%         % 
%         % % Replace image in PSF object with masked image
%         % psfInit.image = masked_psf_image;
% 
%         % If not masking, use this
%         % Generate throwaway PSF object, replace the image with our image
%         par.position = Length([0 0 0], 'nm');
%         psfInit = PSF(par);
%         psfInit.image = psf_image;
% 
%         % Run the fit - given nothing, fmincon - reparameterised
% 
%         % Give it a background noise estimate (let's say you know it to within 10%)
%         if psfInit.backgroundNoise > 1e-5
%             parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
%         else
%             parEst.noiseEstimate = 1e-5;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
%         end
% 
%         parEst.parameterBounds.x = Length([-image_size_nm/2 image_size_nm/2], 'nm');
%         parEst.parameterBounds.y = Length([-image_size_nm/2 image_size_nm/2], 'nm');
%         parEst.parameterBounds.defocus = Length([-10 10], 'nm');
%         parEst.parameterBounds.newangle1 = [-1, 1];
%         parEst.parameterBounds.newangle2 = [-1, 1];
%         parEst.parameterBounds.newangle3 = [0, 1];
% 
%         % Give centroid as initial position estimate
%         parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
%         parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
%         parEst.parameterStartValues.defocus = Length(0, 'nm');
%         parEst.parameterStartValues.newangle1 = 2*rand()-0.5;
%         parEst.parameterStartValues.newangle2 = 2*(rand()-0.5);
%         parEst.parameterStartValues.newangle3 = rand();
% 
%         % Now optimising over photon count too
%         photon_number = par.nPhotons;
%         photon_estimate = round(sum(sum(psfInit.image - parEst.noiseEstimate)));
%         parEst.parameterStartValues.photons = photon_estimate;
%         parEst.parameterBounds.photons = [1, 1e10];%[photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
%         % disp(parEst.parameterBounds.photons);
%         % disp(parEst.parameterStartValues.photons);
%         clear par.nPhotons; % Don't let it know number of photons in advance
% 
%         fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);
% 
%         angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
%         angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
%         angleInclination_estimate = mod(angleInclination_estimate, pi/2);
%         angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);
% 
%         positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
%         positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);
%         photons_estimate = fitResult.estimatesPositionDefocus.ML(7);
%         % disp(photons_estimate)
%         % objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);
% 
%         % % Now we need the value of the objective function at ground truth
%         % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
%         % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
%         % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');
% 
%         % wrap these in abs if doing that again
%         positionX_nm_error = positionX_nm_array(blob_index) - positionX_nm_estimate;
%         positionY_nm_error = positionY_nm_array(blob_index) - positionY_nm_estimate;
%         angleInclination_error = angleInclination_array(blob_index) - angleInclination_estimate;
%         angleAzimuth_error = angleAzimuth_array(blob_index) - angleAzimuth_estimate;
%         photons_error = photon_number - photons_estimate;
%         % objective_function_error = objective_function_true - objective_function_estimate;
% 
%         % Store results
%         positionX_nm_estimate_frame(blob_index) = positionX_nm_estimate;
%         positionY_nm_estimate_frame(blob_index) = positionY_nm_estimate;
%         angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
%         angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
%         photons_estimates_frame(blob_index) = photons_estimate;
%         positionX_nm_errors_frame(blob_index) = positionX_nm_error;
%         positionY_nm_errors_frame(blob_index) = positionY_nm_error;
%         angleInclination_errors_frame(blob_index) = angleInclination_error;
%         angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
%         photons_errors_frame(blob_index) = photons_error;
%         % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
%         % objective_function_estimates_frame(blob_index) = 999999;%objective_function_estimate;
%         % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;
% 
%     end % end loop over blobs
% 
%     % % Show patches on image
%     % figure;
%     % imshow(outlined_image);
%     % title('Original Image with Patch Outlines');
% 
%     positionX_nm_true{frame_index} = positionX_nm_array;
%     positionY_nm_true{frame_index} = positionY_nm_array;
%     angleInclination_true{frame_index} = angleInclination_array;
%     angleAzimuth_true{frame_index} = angleAzimuth_array;
%     photons_true{frame_index} = photon_number * ones(size(angleAzimuth_array));
% 
%     positionX_nm_estimates{frame_index} = positionX_nm_estimate_frame;
%     positionY_nm_estimates{frame_index} = positionY_nm_estimate_frame;
%     angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
%     angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
%     photons_estimates{frame_index} = photons_estimates_frame;
% 
%     positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
%     positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
%     angleInclination_errors{frame_index} = angleInclination_errors_frame;
%     angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
%     photons_errors{frame_index} = photons_errors_frame;
% 
%     % objective_function_trues{frame_index} = objective_function_trues_frame;
%     % objective_function_estimates{frame_index} = objective_function_estimates_frame;
%     % objective_function_errors{frame_index} = objective_function_errors_frame;
% 
% 
%     frame_time = toc;
%     fprintf('  Time for this frame: %.2f seconds\n', frame_time);
% 
%     % Output
% 
%     % Flatten the array (convert cell array to a numeric array)
%     positionX_nm_true_flat = [positionX_nm_true{:}];
%     positionY_nm_true_flat = [positionY_nm_true{:}];
%     angleInclination_true_flat = [angleInclination_true{:}];
%     angleAzimuth_true_flat = [angleAzimuth_true{:}];
%     photons_true_flat = [photons_true{:}];
% 
%     positionX_nm_estimates_flat = [positionX_nm_estimates{:}];
%     positionY_nm_estimates_flat = [positionY_nm_estimates{:}];
%     angleInclination_estimates_flat = [angleInclination_estimates{:}];
%     angleAzimuth_estimates_flat = [angleAzimuth_estimates{:}];
%     photons_estimates_flat = [photons_estimates{:}];
% 
%     positionX_nm_errors_flat = [positionX_nm_errors{:}];
%     positionY_nm_errors_flat = [positionY_nm_errors{:}];
%     angleInclination_errors_flat = [angleInclination_errors{:}];
%     angleAzimuth_errors_flat = [angleAzimuth_errors{:}];
%     photons_errors_flat = [photons_errors{:}];
% 
%     % objective_function_trues_flat = [objective_function_trues{:}];
%     % objective_function_estimates_flat = [objective_function_estimates{:}];
%     % objective_function_errors_flat = [objective_function_errors{:}];
% 
%     % Save data about this stack of images (all frames compressed into one long array)
% 
%     fileID = fopen(fitting_results_path, 'w');
%     fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.4f, ', positionX_nm_true_flat(1:end)));
%     fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.4f, ', positionY_nm_true_flat(1:end)));
%     fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.4f, ', angleInclination_true_flat(1:end)));
%     fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.4f, ', angleAzimuth_true_flat(1:end)));
%     fprintf(fileID, 'x_est = [%s]\n', sprintf('%.4f, ', positionX_nm_estimates_flat(1:end)));
%     fprintf(fileID, 'y_est = [%s]\n', sprintf('%.4f, ', positionY_nm_estimates_flat(1:end)));
%     fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.4f, ', angleInclination_estimates_flat(1:end)));
%     fprintf(fileID, 'az_est = [%s]\n', sprintf('%.4f, ', angleAzimuth_estimates_flat(1:end)));
%     fprintf(fileID, 'x_err = [%s]\n', sprintf('%.4f, ', positionX_nm_errors_flat(1:end)));
%     fprintf(fileID, 'y_err = [%s]\n', sprintf('%.4f, ', positionY_nm_errors_flat(1:end)));
%     fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.4f, ', angleInclination_errors_flat(1:end)));
%     fprintf(fileID, 'az_err = [%s]\n', sprintf('%.4f, ', angleAzimuth_errors_flat(1:end)));
%     % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
%     % fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
%     % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
%     fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_true_flat(1:end)));
%     fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
%     fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
%     fclose(fileID);
% 
% end % end loop over frames
% 













%% Patch version


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
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/sims_9spot_stitched/';

% Output file
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/results_9spot_stitched.py';

% Model to use
model = 'hinterer';



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




% Locate all the images and settings files in the dir
files = dir(frames_dir);
valid_extensions_image = {'.tif'};%{'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
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

% objective_function_trues = [];
% objective_function_estimates = [];
% objective_function_errors = [];

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
    patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

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
    psf_image = imread(frame_path); % if tif
    psf_image = double(psf_image); % Convert to double for calculations if tif
    % psf_image = readmatrix(frame_path); % if csv

    % % Clip values just if want to display
    % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
    % imshow(display_image)

    % psf_image = psf_image ./ norm(psf_image);

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

        patch_centre_x_nm = positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
        patch_centre_y_nm = positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

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

        % If cropping, use this

        % Use ground truth + artificial error as patch centre
        patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_size_nm, image_size_px, 'x');
        patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_size_nm, image_size_px, 'y');

        % Calculate patch size in pixels (ensure odd number for centered patch)
        patch_width_px = floor(patch_width_nm / pixel_size_nm);
        if mod(patch_width_px, 2) == 0
            patch_width_px = patch_width_px + 1; % Ensure odd number
        end

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
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_size_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_size_px, 'y');

        % Calculate expected position relative to patch center - beause of
        % the boundary checking stuff, it might not quite be zero anymore
        expected_offset_x = positionX_nm_array(blob_index) - actual_patch_center_x_nm;
        expected_offset_y = positionY_nm_array(blob_index) - actual_patch_center_y_nm;

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
            parEst.noiseEstimate = psfInit.backgroundNoise;% * (1 + 0.2 * (rand() - 0.5)); % that additional bit to simulate not knowing it exactly
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
        % disp(parEst.parameterBounds.photons);
        % disp(parEst.parameterStartValues.photons);
        clear par.nPhotons; % Don't let it know number of photons in advance

        fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

        angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
        angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
        angleInclination_estimate = mod(angleInclination_estimate, pi/2);
        angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position
        photons_estimate = fitResult.estimatesPositionDefocus.ML(7);
        % disp(photons_estimate)
        % objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

        % % Now we need the value of the objective function at ground truth
        % ground_truth_params = [positionX_nm_array(blob_index), positionY_nm_array(blob_index), 0, angleInclination_array(blob_index), angleAzimuth_array(blob_index)];
        % ground_truth_PSF = createFitPSF(psfInit, ground_truth_params); 
        % objective_function_true = -sum(psf_image .* log(ground_truth_PSF) - ground_truth_PSF - log(gamma(psf_image + 1)), 'all');

        % wrap these in abs if doing that again
        positionX_nm_error = positionX_nm_array(blob_index) - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_array(blob_index) - positionY_nm_estimate;
        angleInclination_error = angleInclination_array(blob_index) - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_array(blob_index) - angleAzimuth_estimate;
        photons_error = photon_number - photons_estimate;
        % objective_function_error = objective_function_true - objective_function_estimate;

        % Store results
        positionX_nm_estimate_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimate_frame(blob_index) = positionY_nm_estimate;
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        photons_estimates_frame(blob_index) = photons_estimate;
        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
        photons_errors_frame(blob_index) = photons_error;
        % objective_function_trues_frame(blob_index) = 999999;%objective_function_true;
        % objective_function_estimates_frame(blob_index) = 999999;%objective_function_estimate;
        % objective_function_errors_frame(blob_index) = 999999;%objective_function_error;





        fprintf('\n--- Coordinate Tracking (Spot %d) ---\n', blob_index);
        fprintf('Ground truth (nm): (%.2f, %.2f)\n', positionX_nm_array(blob_index), positionY_nm_array(blob_index));

        % Convert ground truth to pixels
        gt_x_px = nm_to_px(positionX_nm_array(blob_index), pixel_size_nm, image_size_px, 'x');
        gt_y_px = nm_to_px(positionY_nm_array(blob_index), pixel_size_nm, image_size_px, 'y');
        fprintf('Ground truth (px): (%.2f, %.2f)\n', gt_x_px, gt_y_px);

        % Show patch info
        fprintf('Patch boundaries (px): (%.2f, %.2f) to (%.2f, %.2f)\n', ...
                patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        fprintf('Patch dimensions (px): %dx%d\n', actual_patch_width_px, actual_patch_height_px);

        % Show actual patch center
        fprintf('Patch center (px): (%.2f, %.2f)\n', actual_patch_center_x_px, actual_patch_center_y_px);
        fprintf('Patch center (nm): (%.2f, %.2f)\n', actual_patch_center_x_nm, actual_patch_center_y_nm);

        % Calculate ground truth relative to patch center (this should be the target for fitting)
        fprintf('Ground truth relative to patch center (nm): (%.2f, %.2f)\n', ...
                positionX_nm_array(blob_index) - actual_patch_center_x_nm, ...
                positionY_nm_array(blob_index) - actual_patch_center_y_nm);

        % After fitting, show the results
        fprintf('Fit result in patch coords (nm): (%.2f, %.2f)\n', ...
                fitResult.estimatesPositionDefocus.ML(1), fitResult.estimatesPositionDefocus.ML(2));
        fprintf('Final estimate in global coords (nm): (%.2f, %.2f)\n', ...
                positionX_nm_estimate, positionY_nm_estimate);
        fprintf('Error (nm): (%.2f, %.2f)\n', ...
                positionX_nm_error, positionY_nm_error);
        fprintf('Error (px): (%.2f, %.2f)\n', ...
                positionX_nm_error/pixel_size_nm, positionY_nm_error/pixel_size_nm);





    end % end loop over blobs

    % % Show patches on image
    % figure;
    % imshow(outlined_image);
    % title('Original Image with Patch Outlines');

    positionX_nm_true{frame_index} = positionX_nm_array;
    positionY_nm_true{frame_index} = positionY_nm_array;
    angleInclination_true{frame_index} = angleInclination_array;
    angleAzimuth_true{frame_index} = angleAzimuth_array;
    photons_true{frame_index} = photon_number * ones(size(angleAzimuth_array));

    positionX_nm_estimates{frame_index} = positionX_nm_estimate_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimate_frame;
    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    photons_estimates{frame_index} = photons_estimates_frame;

    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;
    photons_errors{frame_index} = photons_errors_frame;

    % objective_function_trues{frame_index} = objective_function_trues_frame;
    % objective_function_estimates{frame_index} = objective_function_estimates_frame;
    % objective_function_errors{frame_index} = objective_function_errors_frame;


    frame_time = toc;
    fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    % Output

    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_true_flat = [positionX_nm_true{:}];
    positionY_nm_true_flat = [positionY_nm_true{:}];
    angleInclination_true_flat = [angleInclination_true{:}];
    angleAzimuth_true_flat = [angleAzimuth_true{:}];
    photons_true_flat = [photons_true{:}];

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
    % objective_function_estimates_flat = [objective_function_estimates{:}];
    % objective_function_errors_flat = [objective_function_errors{:}];

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
    % fprintf(fileID, 'obj_tru = [%s]\n', sprintf('%.4f, ', objective_function_trues_flat(1:end)));
    % fprintf(fileID, 'obj_est = [%s]\n', sprintf('%.4f, ', objective_function_estimates_flat(1:end)));
    % fprintf(fileID, 'obj_err = [%s]\n', sprintf('%.4f, ', objective_function_errors_flat(1:end)));
    fprintf(fileID, 'photon_tru = [%s]\n', sprintf('%.4f, ', photons_true_flat(1:end)));
    fprintf(fileID, 'photon_est = [%s]\n', sprintf('%.4f, ', photons_estimates_flat(1:end)));
    fprintf(fileID, 'photon_err = [%s]\n', sprintf('%.4f, ', photons_errors_flat(1:end)));
    fclose(fileID);

end % end loop over frames
