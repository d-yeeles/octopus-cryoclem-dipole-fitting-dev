function fit_stack_commandline_thunderstorm(input_image_path, input_thunderstorm_path, output_results_path, model, patch_width_nm, first_frame, last_frame)

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
    if nargin < 5
        error('Not enough input arguments. Usage: fit_multiple_psfs(input_image_and_thunderstorm_path, model, patch_width_nm, starting_frame_index, ending_frame_index)');
    end

    stack_path = input_image_path;
    thunderstorm_results_path = input_thunderstorm_path;
    % thunderstorm_protocol_path = [input_image_and_thunderstorm_path 'thunderstorm_protocol.csv'];
    results_path = output_results_path;

    % Pull out info about image stack
    stack_info = imfinfo(stack_path);
    image_width_px = stack_info(1).Width;
    image_height_px = stack_info(1).Height;
    number_of_frames = numel(stack_info);

    % !!! read in the data rather than hard-code !!!

    % % Read in experimental params from thunderstorm protocol file
    % if ~exist(thunderstorm_protocol_path, 'file')
    %     error('Could not find thunderstorm protocol CSV file: %s', params_path);
    % end
    % fprintf('Reading parameters from: %s\n', thunderstorm_protocol_path);

    % Experimental parameters
    pixel_size_nm = 52;
    par = struct();
    par.wavelength = Length(500, 'nm'); 
    par.pixelSize = Length(pixel_size_nm, 'nm');
    par.objectiveNA = 2.17;
    par.objectiveFocalLength = Length(770, 'mu');
    par.refractiveIndices = [1.31, 2.17, 2.17];
    par.nDiscretizationBFP = 129;
    par.pixelSensitivityMask = PixelSensitivity.uniform(9);

    % Read in thunderstorm localisations
    thunderstorm_results = readtable(thunderstorm_results_path, 'VariableNamingRule', 'preserve');

    % Extract the required columns
    % Check if optional columns exist - only there in sims
    if ismember("angleInclination", thunderstorm_results.Properties.VariableNames)
        thunderstorm_frame_array = thunderstorm_results.("frame_index");
        thunderstorm_x_array = thunderstorm_results.("dipole_posX_nm");
        thunderstorm_y_array = thunderstorm_results.("dipole_posY_nm"); 
        thunderstorm_inc_array = thunderstorm_results.("angleInclination");
        thunderstorm_az_array = thunderstorm_results.("angleAzimuth");
    else
        thunderstorm_frame_array = thunderstorm_results.("frame");
        thunderstorm_x_array = thunderstorm_results.("x [nm]");
        thunderstorm_y_array = thunderstorm_results.("y [nm]"); 
    end
    

    image_width_nm = image_width_px*pixel_size_nm;
    image_height_nm = image_height_px*pixel_size_nm;

    % % ROI location in nm, with origin at centre of image
    % % (use whole image if ROi width or height = 0)
    % if ROI_width_nm == 0 || ROI_height_nm == 0
    %     ROI_min_x_nm = -image_width_nm/2;
    %     ROI_max_x_nm = image_width_nm/2;
    %     ROI_min_y_nm = -image_height_nm/2;
    %     ROI_max_y_nm = image_height_nm/2;
    % else
    %     ROI_min_x_nm = ROI_centre_x_nm - ROI_width_nm/2;
    %     ROI_max_x_nm = ROI_centre_x_nm + ROI_width_nm/2;
    %     ROI_min_y_nm = ROI_centre_y_nm - ROI_height_nm/2;
    %     ROI_max_y_nm = ROI_centre_y_nm + ROI_height_nm/2;
    % end

    % Loop over each frame
    for frame_index = first_frame:last_frame
        
        tic; % timing each frame
    
        fprintf('----------\n');
        fprintf('FRAME %d/%d\n', frame_index, number_of_frames);
        
        % Load frame
        psf_image = imread(stack_path, frame_index);
        psf_image = double(psf_image); % Convert to double for calculations

        par.backgroundNoise = 0; % !!! what should this be !!!
        par.nPhotons = 2000; % !!! what should this be !!!

        % Take thunderstorm results for this frame
        % Note: thunderstorm origin is top-left corner, so need to adjust to centre
        % these will only be there in sims
        % use ground truth only if sims, otherwise use thunderstorm
        current_frame_mask = thunderstorm_frame_array == frame_index; % Create logical mask for when thunderstorm frame number is equal to current frame index
        if exist('thunderstorm_inc_array', 'var') && exist('thunderstorm_az_array', 'var')
            current_frame_x_array = thunderstorm_x_array(current_frame_mask);
            current_frame_y_array = thunderstorm_y_array(current_frame_mask);
            current_frame_inc_array = thunderstorm_inc_array(current_frame_mask);
            current_frame_az_array = thunderstorm_az_array(current_frame_mask);
        else
            current_frame_x_array = thunderstorm_x_array(current_frame_mask) - image_width_nm/2;
            current_frame_y_array = -(thunderstorm_y_array(current_frame_mask) - image_height_nm/2);
        end

        % % Filter to keep only those within ROI
        % in_roi_mask = current_frame_x_array >= ROI_min_x_nm & ...
        %           current_frame_x_array <= ROI_max_x_nm & ...
        %           current_frame_y_array >= ROI_min_y_nm & ...
        %           current_frame_y_array <= ROI_max_y_nm;
        % current_frame_x_array = current_frame_x_array(in_roi_mask);
        % current_frame_y_array = current_frame_y_array(in_roi_mask);

        % % Clip values just if want to display
        % display_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:)));
        % imshow(display_image)

        % Loop over each blob in a frame,
        % cropping out each blob
        % and running the fit on the cropped image

        for blob_index = 1:length(current_frame_x_array)

            fprintf(' ∘ Blob %d/%d\n', blob_index, length(current_frame_x_array));

            patch_width_px = patch_width_nm/pixel_size_nm;
            patch_width_px = 2*floor(patch_width_px/2) + 1; % This enforces odd number of pixels, as required by Hinterer

            par.nPixels = patch_width_px;

            % Initialise results file if first dipole
            if frame_index == 1 && blob_index == 1
                fileID = fopen(results_path, 'w');
                fprintf(fileID, ['frame_index,dipole_index,x_tru,y_tru,newangle1_tru,newangle2_tru,newangle3_tru,inc_tru,az_tru,x_est,y_est,newangle1_est,newangle2_est,newangle3_est,inc_est,az_est,' ...
                                 'x_err,y_err,newangle1_err,newangle2_err,newangle3_err,inc_err,az_err,photon_tru,photon_est,photon_err,dotProduct_err,obj_est,' ...
                                 'covariance_spher,covariance_cart\n']);
                fclose(fileID);
                fprintf('Created new CSV file: %s\n', results_path);
            end

            patch_centre_x_nm = current_frame_x_array(blob_index);
            patch_centre_y_nm = current_frame_y_array(blob_index);

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
                continue
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

            % Run the fit

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
            parEst.parameterBounds.photons = [photon_estimate*0.5, photon_estimate*1.5];%photon_number + photon_number*0.1]; % true value +/- 10%
            clear par.nPhotons; % Don't let it know number of photons in advance

            fitResult = FitPSF_ML_reparam2(psfInit, parEst, model);

            if strcmpi(model, 'gaussian')

                angleInclination_estimate = 0;
                angleAzimuth_estimate = 0;
                photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(4);
                objective_function_estimate = fitResult.estimatesPositionDefocus.ML(5);

            else

                newangle1_estimate = fitResult.estimatesPositionDefocus.ML(4);
                newangle2_estimate = fitResult.estimatesPositionDefocus.ML(5);
                newangle3_estimate = fitResult.estimatesPositionDefocus.ML(6);
                angleInclination_estimate = acos(newangle3_estimate);
                angleAzimuth_estimate = atan2(newangle2_estimate, newangle1_estimate);
                angleInclination_estimate = mod(angleInclination_estimate, pi/2);
                angleAzimuth_estimate = mod(angleAzimuth_estimate, 2*pi);
                photons_fit_estimate = fitResult.estimatesPositionDefocus.ML(7);
                objective_function_estimate = fitResult.estimatesPositionDefocus.ML(8);

            end


            positionX_nm_estimate_local = fitResult.estimatesPositionDefocus.ML(1);
            positionX_nm_estimate = positionX_nm_estimate_local + actual_patch_center_x_nm; % Convert back to global position
            positionY_nm_estimate_local = fitResult.estimatesPositionDefocus.ML(2);
            positionY_nm_estimate = positionY_nm_estimate_local + actual_patch_center_y_nm; % Convert back to global position
            defocus_estimate = fitResult.estimatesPositionDefocus.ML(3);

            % Finding errors in the simulated data is a bit more complicated
            % because of things like multiple detections etc.
            % so for each localisation, find nearest ground truth to the
            % current patch and pick the smallest one. good enough for now.

            % Find index of closest blob to patch centre and use that as ground truth
            distances_to_ground_truth = sqrt((current_frame_x_array - actual_patch_center_x_nm).^2 + (current_frame_y_array - actual_patch_center_y_nm).^2);
            [min_distance, min_index] = min(distances_to_ground_truth);

            % Use this nearest blob as ground truth for this localisation
            % use ground truth only if sims, otherwise use nearest
            if exist('current_frame_inc_array', 'var') && exist('current_frame_az_array', 'var')
                positionX_nm_thunderstorm = current_frame_x_array(blob_index);
                positionY_nm_thunderstorm = current_frame_y_array(blob_index);
                angleInclination_thunderstorm = current_frame_inc_array(blob_index);
                angleAzimuth_thunderstorm = current_frame_az_array(blob_index);
                newangle1_thunderstorm = sin(angleInclination_thunderstorm)*cos(angleAzimuth_thunderstorm);
                newangle2_thunderstorm = sin(angleInclination_thunderstorm)*sin(angleAzimuth_thunderstorm);
                newangle3_thunderstorm = cos(angleInclination_thunderstorm);
                photons_thunderstorm = photon_estimate;
            else
                positionX_nm_thunderstorm = current_frame_x_array(min_index);
                positionY_nm_thunderstorm = current_frame_y_array(min_index);
                angleInclination_thunderstorm = 0; % because thunderstorm didn't give an estimate
                angleAzimuth_thunderstorm = 0;
                photons_thunderstorm = photon_estimate;
            end

            % disp(abs(current_frame_x_array(blob_index) - current_frame_x_array(min_index)));

            positionX_nm_diff = positionX_nm_thunderstorm - positionX_nm_estimate;
            positionY_nm_diff = positionY_nm_thunderstorm - positionY_nm_estimate;
            newangle1_diff = newangle1_thunderstorm - newangle1_estimate;
            newangle2_diff = newangle2_thunderstorm - newangle2_estimate;
            newangle3_diff = newangle3_thunderstorm - newangle3_estimate;
            angleInclination_diff = angleInclination_thunderstorm - angleInclination_estimate;
            angleAzimuth_diff = angleAzimuth_thunderstorm - angleAzimuth_estimate;
            photons_diff = photons_thunderstorm - photons_fit_estimate;

            % Use opening angle dot product thing instead
            dot_product = newangle1_thunderstorm * newangle1_estimate + ...
                          newangle2_thunderstorm * newangle2_estimate + ...
                          newangle3_thunderstorm * newangle3_estimate;
            abs_dot_product = abs(dot_product);
            % Ensure the dot product is within valid range for arccos
            % (numerical precision can sometimes cause values slightly outside [-1,1])
            abs_dot_product = min(max(abs_dot_product, 0), 1);
            % Calculate the error angle in radians
            dotProduct_diff = acos(abs_dot_product);


            % Append results for each blob to an array for this frame
            positionX_nm_thunderstorms_frame(blob_index) = positionX_nm_thunderstorm;
            positionY_nm_thunderstorms_frame(blob_index) = positionY_nm_thunderstorm;
            newangle1_thunderstorms_frame(blob_index) = newangle1_thunderstorm;
            newangle2_thunderstorms_frame(blob_index) = newangle2_thunderstorm;
            newangle3_thunderstorms_frame(blob_index) = newangle3_thunderstorm;
            angleInclination_thunderstorms_frame(blob_index) = angleInclination_thunderstorm;
            angleAzimuth_thunderstorms_frame(blob_index) = angleAzimuth_thunderstorm;
            photons_thunderstorms_frame(blob_index) = photons_thunderstorm;

            positionX_nm_estimates_frame(blob_index) = positionX_nm_estimate;
            positionY_nm_estimates_frame(blob_index) = positionY_nm_estimate;
            newangle1_estimates_frame(blob_index) = newangle1_estimate;
            newangle2_estimates_frame(blob_index) = newangle2_estimate;
            newangle3_estimates_frame(blob_index) = newangle3_estimate;
            angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
            angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
            photons_estimates_frame(blob_index) = photons_fit_estimate;

            positionX_nm_diffs_frame(blob_index) = positionX_nm_diff;
            positionY_nm_diffs_frame(blob_index) = positionY_nm_diff;
            newangle1_diffs_frame(blob_index) = newangle1_diff;
            newangle2_diffs_frame(blob_index) = newangle2_diff;
            newangle3_diffs_frame(blob_index) = newangle3_diff;
            angleInclination_diffs_frame(blob_index) = angleInclination_diff;
            angleAzimuth_diffs_frame(blob_index) = angleAzimuth_diff;
            photons_diffs_frame(blob_index) = photons_diff;
            dotProduct_diffs_frame(blob_index) = dotProduct_diff;
            objective_function_estimates_frame(blob_index) = objective_function_estimate;






            % % Covariance matrix
            % try
            % 
            %     [covarMatrix, fisherMatrix] = calculateCovarianceMatrix(psfInit, positionX_nm_estimate_local, positionY_nm_estimate_local, defocus_estimate, angleInclination_estimate, angleAzimuth_estimate, photons_fit_estimate, parEst.noiseEstimate, model);
            % 
            %     % Extract standard deviations (square root of diagonal elements)
            %     param_stds = sqrt(diag(covarMatrix));
            % 
            %     % Get uncertainty values based on model type
            %     if strcmpi(model, 'gaussian')
            %         x_std = param_stds(1);
            %         y_std = param_stds(2);
            %         % z_std = param_stds(3);
            %         % photon_std = param_stds(4);
            %         theta_std = NaN;
            %         phi_std = NaN;
            %     else
            %         x_std = param_stds(1);
            %         y_std = param_stds(2);
            %         % z_std = param_stds(3);
            %         theta_std = param_stds(3);
            %         phi_std = param_stds(4);
            %         % photon_std = param_stds(6);
            %     end
            % 
            % catch ME
            %     % Handle any errors in the covariance calculation gracefully
            %     fprintf('  Warning: Could not calculate uncertainties. Error: %s\n', ME.message);
            %     x_std = NaN;
            %     y_std = NaN;
            %     % z_std = NaN;
            %     theta_std = NaN;
            %     phi_std = NaN;
            %     % photon_std = NaN;
            % end
            
            % Covariance matrix
            % [covarMatrix, fisherMatrix] = calculateCovarianceMatrix(psfInit, positionX_nm_estimate_local, positionY_nm_estimate_local, defocus_estimate, newangle1_estimate, newangle2_estimate, newangle3_estimate, photons_fit_estimate, parEst.noiseEstimate, model);
            % [covarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(psfInit, positionX_nm_estimate_local, positionY_nm_estimate_local, defocus_estimate, angleInclination_estimate, angleAzimuth_estimate, photons_fit_estimate, parEst.noiseEstimate, model);
            [cartesianCovarMatrix, sphericalCovarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(psfInit, positionX_nm_estimate_local, positionY_nm_estimate_local, defocus_estimate, newangle1_estimate, newangle2_estimate, newangle3_estimate, photons_fit_estimate, parEst.noiseEstimate, model);

            % % Extract standard deviations (square root of diagonal elements)
            % param_stds = sqrt(diag(covarMatrix));
            % photon_std = param_stds(5);
            % disp(photon_std)



            % Create the matrix string with standard nested array notation for spherical covariance
            covar_str_spher = '[';  % Start with matrix opening
            
            for i = 1:size(sphericalCovarMatrix, 1)
                covar_str_spher = [covar_str_spher, '['];  % Add row opening bracket
                
                for j = 1:size(sphericalCovarMatrix, 2)
                    % Format the number in scientific notation
                    covar_str_spher = [covar_str_spher, sprintf('%.6e', sphericalCovarMatrix(i,j))];
                    
                    % Add comma between elements unless it's the last element in the row
                    if j < size(sphericalCovarMatrix, 2)
                        covar_str_spher = [covar_str_spher, '|'];
                    end
                end
                
                covar_str_spher = [covar_str_spher, ']'];  % Close this row
                
                % Add comma between rows unless it's the last row
                if i < size(sphericalCovarMatrix, 1)
                    covar_str_spher = [covar_str_spher, '|'];
                end
            end
            
            covar_str_spher = [covar_str_spher, ']'];  % Close the matrix
            


            
            % Create the matrix string with standard nested array notation for cartesian covariance
            covar_str_cart = '[';  % Start with matrix opening
            
            for i = 1:size(cartesianCovarMatrix, 1)
                covar_str_cart = [covar_str_cart, '['];  % Add row opening bracket
                
                for j = 1:size(cartesianCovarMatrix, 2)
                    % Format the number in scientific notation
                    covar_str_cart = [covar_str_cart, sprintf('%.6e', cartesianCovarMatrix(i,j))];
                    
                    % Add comma between elements unless it's the last element in the row
                    if j < size(cartesianCovarMatrix, 2)
                        covar_str_cart = [covar_str_cart, '|'];
                    end
                end
                
                covar_str_cart = [covar_str_cart, ']'];  % Close this row
                
                % Add comma between rows unless it's the last row
                if i < size(cartesianCovarMatrix, 1)
                    covar_str_cart = [covar_str_cart, '|'];
                end
            end
            
            covar_str_cart = [covar_str_cart, ']'];  % Close the matrix




            fileID = fopen(results_path, 'a');
            
            % Format all the numeric values first
            numeric_part = sprintf('%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
                frame_index, ...
                blob_index, ...
                positionX_nm_thunderstorm, ...
                positionY_nm_thunderstorm, ...
                newangle1_thunderstorm, ...
                newangle2_thunderstorm, ...
                newangle3_thunderstorm, ...
                angleInclination_thunderstorm, ...
                angleAzimuth_thunderstorm, ...
                positionX_nm_estimate, ...
                positionY_nm_estimate, ...
                newangle1_estimate, ...
                newangle2_estimate, ...
                newangle3_estimate, ...
                angleInclination_estimate, ...
                angleAzimuth_estimate, ...
                positionX_nm_diff, ...
                positionY_nm_diff, ...
                newangle1_diff, ...
                newangle2_diff, ...
                newangle3_diff, ...
                angleInclination_diff, ...
                angleAzimuth_diff, ...
                photons_thunderstorm, ...
                photons_fit_estimate, ...
                photons_diff, ...
                dotProduct_diff, ...
                objective_function_estimate);
            
            % Write everything to file
            fprintf(fileID, '%s,%s,%s\n', numeric_part, covar_str_spher, covar_str_cart);
            fclose(fileID);


        end % end loop over blobs





        % % -----------------------------------------------------------------
        % % VISUALISATION
        % % -----------------------------------------------------------------
        % % Display the original frame image
        % figure(100); clf;
        % imagesc(psf_image);
        % colormap(gray);
        % axis image;
        % hold on;
        % 
        % % Plot ROI in yellow
        % ROI_min_x_px = nm_to_px(ROI_min_x_nm, pixel_size_nm, image_width_px, 'x');
        % ROI_max_x_px = nm_to_px(ROI_max_x_nm, pixel_size_nm, image_width_px, 'x');
        % ROI_min_y_px = nm_to_px(ROI_min_y_nm, pixel_size_nm, image_width_px, 'y');
        % ROI_max_y_px = nm_to_px(ROI_max_y_nm, pixel_size_nm, image_width_px, 'y');
        % % Make sure left is smaller than right, and top is smaller than bottom
        % roi_left = min(ROI_min_x_px, ROI_max_x_px);
        % roi_top = min(ROI_min_y_px, ROI_max_y_px);
        % roi_width = abs(ROI_max_x_px - ROI_min_x_px);
        % roi_height = abs(ROI_max_y_px - ROI_min_y_px);
        % 
        % % Draw the rectangle using corrected values
        % rectangle('Position', [roi_left, roi_top, roi_width, roi_height], ...
        %     'EdgeColor', 'y', 'LineWidth', 1);
        % % rectangle('Position', [ROI_min_x_px, ROI_min_y_px, ROI_max_x_px-ROI_min_x_px, ROI_max_y_px-ROI_min_y_px], ...
        % %     'EdgeColor', 'y', 'LineWidth', 1);
        % 
        % % Plot all patch boundaries in red
        % for b_idx = 1:length(current_frame_x_array)
        %     patch_centre_x_px = nm_to_px(current_frame_x_array(b_idx), pixel_size_nm, image_width_px, 'x');
        %     patch_centre_y_px = nm_to_px(current_frame_y_array(b_idx), pixel_size_nm, image_height_px, 'y');
        % 
        %     % Half width for calculations
        %     half_width = floor(patch_width_px / 2);
        % 
        %     % Calculate patch boundaries
        %     patch_start_x = max(1, round(patch_centre_x_px - half_width));
        %     patch_start_y = max(1, round(patch_centre_y_px - half_width));
        %     patch_end_x = min(image_width, patch_start_x + patch_width_px - 1);
        %     patch_end_y = min(image_height, patch_start_y + patch_width_px - 1);
        % 
        %     % Draw rectangle around the patch
        %     rectangle('Position', [patch_start_x, patch_start_y, patch_end_x-patch_start_x, patch_end_y-patch_start_y], ...
        %         'EdgeColor', 'r', 'LineWidth', 1);
        % end
        % 
        % % Plot thunderstorm localizations as green dots
        % for ts_idx = 1:length(current_frame_x_array)
        %     ts_x_px = nm_to_px(positionX_nm_thunderstorms_frame(ts_idx), pixel_size_nm, image_width_px, 'x');
        %     ts_y_px = nm_to_px(positionY_nm_thunderstorms_frame(ts_idx), pixel_size_nm, image_height_px, 'y');
        %     plot(ts_x_px, ts_y_px, 'g.', 'MarkerSize', 1);
        % end
        % 
        % % Plot estimated localizations as yellow dots
        % for est_idx = 1:length(current_frame_x_array)
        %     est_x_px = nm_to_px(positionX_nm_estimates_frame(est_idx), pixel_size_nm, image_width_px, 'x');
        %     est_y_px = nm_to_px(positionY_nm_estimates_frame(est_idx), pixel_size_nm, image_height_px, 'y');
        %     plot(est_x_px, est_y_px, 'y.', 'MarkerSize', 1);
        % end
        % 
        % % % Add a legend
        % % legend('', 'Patch boundaries', 'ThunderStorm', 'Estimated', 'Location', 'NorthEast');
        % title(sprintf('Frame %d - Localizations', frame_index));
        % 
        % % Add a colorbar
        % colorbar;
        % 
        % % Force the figure to update
        % drawnow;
        % 
        % % Optionally save the figure
        % saveas(gcf, sprintf('%sframe_%03d_localization.png', input_image_and_thunderstorm_path, frame_index));
        % 
        % % -----------------------------------------------------------------




        frame_time = toc;
        fprintf('  Time for this frame: %.2f seconds\n', frame_time);

    end % end loop over frames

end




