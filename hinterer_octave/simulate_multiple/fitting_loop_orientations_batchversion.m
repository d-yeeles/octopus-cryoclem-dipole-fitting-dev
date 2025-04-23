function fitting_loop_orientations_batchversion(frames_dir, fitting_results_path, model_name)

    % frames_dir: path to directory of simulation frames to analyse
    % fitting_results_path: path and filename of where results are saved
    % model_name: whether to fit 'hinterer' or 'gaussian' or 'correct_polar or 'correct_orientations'



    % addpath(genpath('../'));
    addpath(genpath('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/'));


    % Input params
    % frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_inc90/';
    
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
    
    % Loop over each frame
    for frame_index = 1:length(frame_paths)
    
        fprintf('----------\n');
        fprintf('FRAME %d\n', frame_index);
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
    
        % Convert to greyscale if not
        if size(psf_image, 3) == 3
            psf_image = rgb2gray(psf_image);
        end
    
        psf_image = double(psf_image); % Convert to double for calculations
        psf_image = max(min(psf_image, 255), 0); % Clip values to [0, 255]
        psf_image = uint8(psf_image); % Convert back to uint8
        psf_image = double(psf_image);
        psf_image = (psf_image - min(psf_image(:))) / (max(psf_image(:)) - min(psf_image(:))); % need to do this to normalise it to 0,1
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
    
            centroid_x_nm = positionX_nm_array(blob_index) + (rand - 0.5) * 1000; % this is simulating some other blob-finding result
            centroid_y_nm = positionY_nm_array(blob_index) + (rand - 0.5) * 1000;
    
            centroid_x_px = positionX_px_array(blob_index);
            centroid_y_px = positionY_px_array(blob_index);
    
            patch_start_x_nm = centroid_x_nm - patch_width_nm/2;
            patch_start_y_nm = centroid_y_nm - patch_width_nm/2;
            patch_end_x_nm = centroid_x_nm + patch_width_nm/2;
            patch_end_y_nm = centroid_y_nm + patch_width_nm/2;
    
            patch_start_x_px = round(centroid_x_px - (patch_width_px + 1)/2);
            patch_start_y_px = round(centroid_y_px - (patch_width_px + 1)/2);
            patch_end_x_px = round(centroid_x_px + (patch_width_px + 1)/2);
            patch_end_y_px = round(centroid_y_px + (patch_width_px + 1)/2);
    
            % Generate throwaway PSF object, later replace the image with our masked image
            par.position = Length([0 0 0], 'nm');
            psfInit = PSF(par);
    
            % Create masked image, where everything outside patch = 0
            masked_psf_image = psf_image;
            [image_height, image_width] = size(masked_psf_image);  % Get the image size
            
            mask_x = (1:image_width) >= patch_start_x_px & (1:image_width) <= patch_end_x_px;  % x indices within the patch
            mask_y = (1:image_height) >= patch_start_y_px & (1:image_height) <= patch_end_y_px;  % y indices within the patch
            
            masked_psf_image(~mask_y, :) = 0;
            masked_psf_image(:, ~mask_x) = 0;
    
            % Replace image in PSF object with masked image
            psfInit.image = masked_psf_image;
           
            % Give patch as localisation bounds
            parEst.parameterBounds.x = Length([patch_start_x_nm patch_end_x_nm], 'nm');
            parEst.parameterBounds.y = Length([patch_start_y_nm patch_end_y_nm], 'nm');
            parEst.parameterBounds.defocus = Length([-10 10], 'nm');
            parEst.parameterBounds.inclination = [0 pi/2];
            parEst.parameterBounds.azimuth = [-Inf Inf];
            % Give centroid as initial position estimate
            parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
            parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
            parEst.parameterStartValues.defocus = Length(0, 'nm');
            parEst.parameterStartValues.inclination = 0; % add angle optimiser
            parEst.parameterStartValues.azimuth = 0; % add angle optimiser
    
            % Run the fit - with angle optimising
            if strcmp(model_name, 'hinterer')
                tic;
                fitResult = FitPSF_angles(psfInit, parEst);
                elapsed_time = toc;
                fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
                angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
                angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);
            end
            if strcmp(model_name, 'gaussian')
                tic;
                fitResult = FitPSF_gaussian(psfInit, parEst);
                elapsed_time = toc;
                fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
                angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
                angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);
            end
            if strcmp(model_name, 'correct_polar')
                parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
                tic;
                fitResult = FitPSF_notpolar(psfInit, parEst);
                elapsed_time = toc;
                fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
                angleInclination_estimate = parEst.angleInclinationEstimate;
                angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(4); 
            end
            if strcmp(model_name, 'correct_orientation')
                parEst.angleInclinationEstimate = angleInclination_array(blob_index); % given orientations in this case
                parEst.angleAzimuthEstimate = angleAzimuth_array(blob_index); % given orientations in this case
                tic;
                fitResult = FitPSF(psfInit, parEst);
                elapsed_time = toc;
                fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
                angleInclination_estimate = parEst.angleInclinationEstimate;
                angleAzimuth_estimate = parEst.angleAzimuthEstimate;
            end


    
            positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
            positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);
    
            % wrap these in abs if doing that again
            positionX_nm_error = positionX_nm_array(blob_index) - positionX_nm_estimate;
            positionY_nm_error = positionY_nm_array(blob_index) - positionY_nm_estimate;
            angleInclination_error = angleInclination_array(blob_index) - angleInclination_estimate;
            angleAzimuth_error = angleAzimuth_array(blob_index) - angleAzimuth_estimate;
    
            % Store results
            positionX_nm_estimate_frame(blob_index) = positionX_nm_estimate;
            positionY_nm_estimate_frame(blob_index) = positionY_nm_estimate;
            angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
            angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
            positionX_nm_errors_frame(blob_index) = positionX_nm_error;
            positionY_nm_errors_frame(blob_index) = positionY_nm_error;
            angleInclination_errors_frame(blob_index) = angleInclination_error;
            angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;
    
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
    
    end % end loop over frames
    
    % Flatten the array (convert cell array to a numeric array)
    positionX_nm_true = [positionX_nm_true{:}];
    positionY_nm_true = [positionY_nm_true{:}];
    angleInclination_true = [angleInclination_true{:}];
    angleAzimuth_true = [angleAzimuth_true{:}];
    
    positionX_nm_estimates = [positionX_nm_estimates{:}];
    positionY_nm_estimates = [positionY_nm_estimates{:}];
    angleInclination_estimates = [angleInclination_estimates{:}];
    angleAzimuth_estimates = [angleAzimuth_estimates{:}];
    
    positionX_nm_errors = [positionX_nm_errors{:}];
    positionY_nm_errors = [positionY_nm_errors{:}];
    angleInclination_errors = [angleInclination_errors{:}];
    angleAzimuth_errors = [angleAzimuth_errors{:}];
    
    % Save data about this stack of images (all frames compressed into one long array)
    % fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_inc90/fitting_results_hinterer.py';
    
    fileID = fopen(fitting_results_path, 'w');
    fprintf(fileID, 'x_tru = [%s]\n', sprintf('%.2f, ', positionX_nm_true(1:end)));
    fprintf(fileID, 'y_tru = [%s]\n', sprintf('%.2f, ', positionY_nm_true(1:end)));
    fprintf(fileID, 'inc_tru = [%s]\n', sprintf('%.2f, ', angleInclination_true(1:end)));
    fprintf(fileID, 'az_tru = [%s]\n', sprintf('%.2f, ', angleAzimuth_true(1:end)));
    fprintf(fileID, 'x_est = [%s]\n', sprintf('%.2f, ', positionX_nm_estimates(1:end)));
    fprintf(fileID, 'y_est = [%s]\n', sprintf('%.2f, ', positionY_nm_estimates(1:end)));
    fprintf(fileID, 'inc_est = [%s]\n', sprintf('%.2f, ', angleInclination_estimates(1:end)));
    fprintf(fileID, 'az_est = [%s]\n', sprintf('%.2f, ', angleAzimuth_estimates(1:end)));
    fprintf(fileID, 'x_err = [%s]\n', sprintf('%.2f, ', positionX_nm_errors(1:end)));
    fprintf(fileID, 'y_err = [%s]\n', sprintf('%.2f, ', positionY_nm_errors(1:end)));
    fprintf(fileID, 'inc_err = [%s]\n', sprintf('%.2f, ', angleInclination_errors(1:end)));
    fprintf(fileID, 'az_err = [%s]\n', sprintf('%.2f, ', angleAzimuth_errors(1:end)));
    fclose(fileID);

end
