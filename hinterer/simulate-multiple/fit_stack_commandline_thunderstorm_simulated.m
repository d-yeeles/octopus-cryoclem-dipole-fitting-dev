function fit_stack_commandline_thunderstorm_simulated(stack_path, params_path, thunderstorm_results_path, results_path, model, patch_width_nm, starting_frame_index, ending_frame_index)
    % Function to fit multiple PSFs in a frame
    % 
    % USAGE:
    %   fit_multiple_psfs(stack_dir, stack_name, results_path, model, patch_width_nm)
    %
    % INPUTS:
    %   stack_dir - Directory containing the image stack and parameter files
    %   stack_path - Full path to the image stack file
    %   results_path - Output path for fitting results
    %   model - Model to use ('mortensen', 'hinterer', or 'gaussian')
    %   patch_width_nm - how big of a patch to consider around each dipole
    %   (suggest 988 nm)

    % Add this after the function signature
    if nargin < 8
        error('Not enough input arguments. Usage: fit_multiple_psfs(stack_path, params_path, thunderstorm_results_path, results_path, model, patch_width_nm, starting_frame_index, ending_frame_index)');
    end

    % Pull out info about image stack
    stack_info = imfinfo(stack_path);
    image_width_px = stack_info(1).Width;
    image_height_px = stack_info(1).Height;
    number_of_frames = numel(stack_info);

    % Find and read the CSV parameters file
    if ~exist(params_path, 'file')
        error('Could not find parameters CSV file: %s', params_path);
    end

    % Read the CSV file
    fprintf('Reading parameters from: %s\n', params_path);
    opts = detectImportOptions(params_path);
    opts.VariableNamingRule = 'preserve';
    params_table = readtable(params_path, opts);

    % Get universal params
    pixel_size_nm = params_table.pixel_size_nm(1);

    % Set up parameters for PSF modeling
    par = struct();
    par.wavelength = Length(params_table.wavelength(1), 'nm'); 
    par.pixelSize = Length(pixel_size_nm, 'nm');
    par.objectiveNA = params_table.objectiveNA(1);
    par.objectiveFocalLength = Length(params_table.objectiveFocalLength(1), 'mu');

    % For refractive indices, combine the three columns
    par.refractiveIndices = str2num(params_table.refractiveIndices{1});

    par.nDiscretizationBFP = params_table.nDiscretizationBFP(1);
    par.backgroundNoise = params_table.backgroundNoise(1);
    par.nPhotons = params_table.nPhotons(1); % won't have this in reality
    par.pixelSensitivityMask = PixelSensitivity.uniform(9);

    % Read in thunderstorm localisations
    thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');

    % Extract the required columns
    thunderstorm_frame_array = thunderstorm_results.("frame");
    thunderstorm_x_array = thunderstorm_results.("x [nm]");
    thunderstorm_y_array = thunderstorm_results.("y [nm]");

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % Loop over each image path and process

    % Loop over each frame in stack
    for frame_index = starting_frame_index:ending_frame_index

        tic; % timing each frame

        fprintf('----------\n');
        fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
        
        % Load frame
        psf_image = imread(stack_path, frame_index);
        psf_image = double(psf_image); % Convert to double for calculations

        % Read in true angle (won't have this in reality of course)
        % Find the row that matches current frame and dipole
        this_frame_indices = find(params_table.frame_index == frame_index);
        positionX_nm_array = params_table.dipole_posX_nm(this_frame_indices);
        positionY_nm_array = params_table.dipole_posY_nm(this_frame_indices);
        angleInclination_array = params_table.angleInclination(this_frame_indices);
        angleAzimuth_array = params_table.angleAzimuth(this_frame_indices);

        % Take thunderstorm results for that frame
        % Note: thunderstorm origin is top-left corner, so need to adjust to centre
        current_frame_mask = thunderstorm_frame_array == frame_index;
        current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
        current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);

        % % Clip values just if want to display
        % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
        % imshow(display_image)

        % Loop over each blob in a frame,
        % cropping out each blob
        % and running the fit on the cropped image

        for blob_index = 1:length(current_frame_x_array)

            fprintf(' ∘ Blob %d/%d\n', blob_index, length(current_frame_x_array));

            positionX_nm_true = current_frame_x_array(blob_index);
            positionY_nm_true = current_frame_y_array(blob_index);

            patch_width_px = patch_width_nm/pixel_size_nm;
            patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

            par.nPixels = patch_width_px;

            % Initialise results file if first dipole
            if frame_index == 1 && blob_index == 1
                fileID = fopen(results_path, 'w');
                fprintf(fileID, 'frame_index,dipole_index,x_tru,y_tru,inc_tru,az_tru,x_est,y_est,inc_est,az_est,x_err,y_err,inc_err,az_err,photon_tru,photon_est,photon_err,obj_est\n');
                fclose(fileID);
                fprintf('Created new CSV file: %s\n', results_path);
            end

            patch_centre_x_nm = positionX_nm_true;%positionX_nm_array(blob_index);% + (rand - 0.5) * 200; % this is simulating some other blob-finding result
            patch_centre_y_nm = positionY_nm_true;%positionY_nm_array(blob_index);% + (rand - 0.5) * 200;

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

            % First calculate the integer patch boundaries that ensure the patch stays within image
            % patch_start_x_px = max(1, floor(patch_centre_x_px - half_width));
            % patch_start_y_px = max(1, floor(patch_centre_y_px - half_width));
            patch_start_x_px = max(1, round(patch_centre_x_px - half_width));
            patch_start_y_px = max(1, round(patch_centre_y_px - half_width));
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

            % Skip if nothing in the patch
            if ~any(patch_psf_image)
                disp('Patch contains nothing. Skipping it.');
                return
            end

            % Calculate expected position relative to patch center
            expected_offset_x = patch_centre_x_nm - actual_patch_center_x_nm;
            expected_offset_y = patch_centre_y_nm - actual_patch_center_y_nm;

            % Setup PSF object for fitting
            par.position = Length([0 0 0], 'nm');
            % psfInit = PSF(par);
            if strcmpi(model, 'hinterer')
                psfInit = PSF(par);
            elseif strcmpi(model, 'mortensen')
                psfInit = PSF_mortensen(par);
            elseif strcmpi(model, 'gaussian')
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
            parEst.parameterBounds.newangle1 = [-1, 1];
            parEst.parameterBounds.newangle2 = [-1, 1];
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

            if strcmpi(model, 'gaussian')

                angleInclination_estimate = 0;
                angleAzimuth_estimate = 0;
                photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(4);
                objective_function_estimate = fitResult.estimatesPositionDefocus.ML(5);

            else

                angleInclination_estimate = acos(fitResult.estimatesPositionDefocus.ML(6));
                angleAzimuth_estimate = atan2(fitResult.estimatesPositionDefocus.ML(5), fitResult.estimatesPositionDefocus.ML(4));
                angleInclination_estimate = mod(angleInclination_estimate, pi/2);
                angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);
                photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
                objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

            end

            positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1) + actual_patch_center_x_nm; % Convert back to global position
            positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2) + actual_patch_center_y_nm; % Convert back to global position

            % Finding errors in the simulated data is a bit more complicated
            % because of things like multiple detections etc.
            % so for each localisation, find nearest ground truth to the
            % current patch and pick the smallest one. good enough for now.

            % Find index of closest blob to patch centre and use that as ground truth
            distances_to_ground_truth = sqrt((positionX_nm_array - actual_patch_center_x_nm).^2 + (positionY_nm_array - actual_patch_center_y_nm).^2);
            [min_distance, min_index] = min(distances_to_ground_truth);

            % Use this nearest blob as ground truth for this localisation
            positionX_nm_true = positionX_nm_array(min_index); % won't have this in reality
            positionY_nm_true = positionY_nm_array(min_index); % won't have this in reality
            angleInclination_true = angleInclination_array(min_index); % won't have this in reality
            angleAzimuth_true = angleAzimuth_array(min_index); % won't have this in reality
            photons_true = photon_estimate;

            positionX_nm_error = positionX_nm_true - positionX_nm_estimate;
            positionY_nm_error = positionY_nm_true - positionY_nm_estimate;
            angleInclination_error = angleInclination_true - angleInclination_estimate;
            angleAzimuth_error = angleAzimuth_true - angleAzimuth_estimate;
            photons_error = photons_true - photons_fit_estimate;

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
            objective_function_estimates_frame(blob_index) = objective_function_estimate;

            % Write the result for this dipole to the CSV file
            fileID = fopen(results_path, 'a');  % Open in append mode
            fprintf(fileID, '%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
                frame_index, ...
                blob_index, ...
                positionX_nm_true, ...
                positionY_nm_true, ...
                angleInclination_true, ...
                angleAzimuth_true, ...
                positionX_nm_estimate, ...
                positionY_nm_estimate, ...
                angleInclination_estimate, ...
                angleAzimuth_estimate, ...
                positionX_nm_error, ...
                positionY_nm_error, ...
                angleInclination_error, ...
                angleAzimuth_error, ...
                photons_true, ...
                photons_fit_estimate, ...
                photons_error, ...
                objective_function_estimate);
            fclose(fileID);

        end % end loop over blobs

        frame_time = toc;
        fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    end % end loop over frames

end




