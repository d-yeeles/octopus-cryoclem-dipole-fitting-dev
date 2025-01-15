%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

% Global params - these will be the same whether sim or fit

output_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output.tif';
number_of_spots = 4;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = sqrt(number_of_spots)*2000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd

% Attocube params
par.nPixels = image_size_px;
par.wavelength = Length(500,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(770,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
par.backgroundNoise = 1/number_of_spots; % taken from looking at blank bit of example data
par.nPhotons = 500;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Initialise empty position arrays for fitting later
% (PSF object wants it in nm with origin at centre, drawing stuff wants px with origin at top-left)
positionX_nm_array = [];
positionY_nm_array = [];
positionX_px_array = [];
positionY_px_array = [];
angleInclination_array = [];
angleAzimuth_array = [];

% Loop over a bunch of random orientations, positions, photon counts
for i = 1:number_of_spots

    % regular grid
    % if number_of_spots == 1
    %     positionX_nm = -1249.4;
    %     positionY_nm = -671.1673;
    % else
    %     % test grid
    %     % uses relative (x,y) where origin is bottom-left and they span 0-1
    %     relative_x = inner_bound + outer_bound * ((mod(i-1, sqrt(number_of_spots))) / (sqrt(number_of_spots) - 1));
    %     relative_y = inner_bound + outer_bound * ((sqrt(number_of_spots) - 1 - floor((i-1) / sqrt(number_of_spots))) / (sqrt(number_of_spots) - 1));
    %     positionX_nm = (relative_x - 0.5) * image_size_nm;
    %     positionY_nm = (relative_y - 0.5) * image_size_nm;
    % end

    % random but avoiding edges and each other

    min_distance_nm = 1000;
    valid_position = false;
    timeout = 20;  % timeout
    start_time = tic;

    while ~valid_position

        % Check if the timeout has been exceeded
        elapsed_time = toc(start_time);
        if elapsed_time > timeout
            disp("Timeout. Probably can't fit that many blobs in this small a space.");
            return;  % Exit the loop if 30 seconds have passed
        end

        % Generate a random position avoiding edges
        relative_x = inner_bound + outer_bound * rand();
        relative_y = inner_bound + outer_bound * rand();

        % Convert to nm position
        positionX_nm = (relative_x - 0.5) * image_size_nm;
        positionY_nm = (relative_y - 0.5) * image_size_nm;

        % Check distance from all existing spots
        if isempty(positionX_nm_array)
            valid_position = true; % First spot is always valid
        else
            distances = sqrt((positionX_nm_array - positionX_nm).^2 + ...
                             (positionY_nm_array - positionY_nm).^2);
            if all(distances >= min_distance_nm)
                valid_position = true;
            end
        end
    end

    % we'll need its pixel equivalent for drawing later on
    positionX_px = nm_to_px(positionX_nm, pixel_size_nm, image_size_px, 'x');
    positionY_px = nm_to_px(positionY_nm, pixel_size_nm, image_size_px, 'y');

    angleInclination = pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
    angleAzimuth = 2*pi*rand;%0; % doesn't affect anything

    par.position = Length([positionX_nm positionY_nm 0], 'nm'); % This is nm away from centre of image
    par.dipole = Dipole(angleInclination, angleAzimuth);

    psf = PSF(par);

    % if first iteration, use this psf image
    % later loops just add to this image
    if i == 1
        psf_total_image = psf.image;
    else
        psf_total_image = psf_total_image + psf.image;
    end

    positionX_nm_array(end+1) = positionX_nm;
    positionY_nm_array(end+1) = positionY_nm;
    positionX_px_array(end+1) = positionX_px;
    positionY_px_array(end+1) = positionY_px;
    angleInclination_array(end+1) = angleInclination;
    angleAzimuth_array(end+1) = angleAzimuth;

end

disp('ground truth')
fprintf('x = [%s]\n', sprintf('%.2f, ', positionX_nm_array(1:end)));
fprintf('y = [%s]\n', sprintf('%.2f, ', positionY_nm_array(1:end)));
fprintf('inc = [%s]\n', sprintf('%.2f, ', angleInclination_array(1:end)*180/pi));
fprintf('az = [%s]\n', sprintf('%.2f, ', angleAzimuth_array(1:end)*180/pi));


% Plot
figure;
subplot(1,3,1);
psf_total_image = mat2gray(psf_total_image);
psf_total_image = uint8(255*psf_total_image);
imshow(psf_total_image)
title('Simulated data');

% Output as tif stack
imwrite(psf_total_image, output_path);
fprintf('Simulation output to \n %s\n', output_path);



%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/';
patch_width_nm = 1000; % size of patch around blob to consider
patch_width_px = patch_width_nm/pixel_size_nm;

% Create the centroids in image coordinates
% 'Read in' ground truth positions (this will be rpelaced by like reading in thunderstorm estimate or something)

f_array = ones(1, number_of_spots);
x_array = positionX_px_array; % replace this with output of some blob-finding algorithm
y_array = positionY_px_array;

% Loop over each frame:
files = dir(frames_dir);
valid_extensions = {'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif'};
frame_paths = {};

% Keep only image files
for i = 1:length(files)
    [~, ~, ext] = fileparts(files(i).name);
    if ismember(lower(ext), valid_extensions)
        frame_paths{end+1} = fullfile(frames_dir, files(i).name);
    end
end

% Loop over each image path and process

% Initialise results arrays
angleInclination_estimates = [];
angleAzimuth_estimates = [];
positionX_estimates = [];
positionY_estimates = [];
defocus_estimates = [];

for frame_index = 1:length(frame_paths)

    fprintf('----------\n');
    fprintf('FRAME %d\n', frame_index);
    frame_path = frame_paths{frame_index};

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
    angleInclination_estimates_frame = [];
    angleAzimuth_estimates_frame = [];
    positionX_nm_estimate_frame = [];
    positionY_nm_estimate_frame = [];
    defocus_estimates = [];
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
        parEst.parameterBounds.inclination = [0 pi];
        parEst.parameterBounds.azimuth = [0 2*pi];
        % Give centroid as initial position estimate
        parEst.parameterStartValues.x = Length(centroid_x_nm, 'nm');
        parEst.parameterStartValues.y = Length(centroid_y_nm, 'nm');
        parEst.parameterStartValues.defocus = Length(0, 'nm');
        parEst.parameterStartValues.inclination = 0; % add angle optimiser
        parEst.parameterStartValues.azimuth = 0; % add angle optimiser

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
        tic;
        fitResult = FitPSF_angles(psfInit, parEst);
        elapsed_time = toc;
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

        positionX_nm_estimate = fitResult.estimatesPositionDefocus.ML(1);
        positionY_nm_estimate = fitResult.estimatesPositionDefocus.ML(2);
        defocus_estimate = fitResult.estimatesPositionDefocus.ML(3);


        % disp(angleInclination_array(blob_index))
        % disp(angleInclination_estimate)

        % wrap these in abs if doing that again
        positionX_nm_error = positionX_nm_array(blob_index) - positionX_nm_estimate;
        positionY_nm_error = positionY_nm_array(blob_index) - positionY_nm_estimate;
        angleInclination_error = angleInclination_array(blob_index) - angleInclination_estimate;
        angleAzimuth_error = angleAzimuth_array(blob_index) - angleAzimuth_estimate;

        % Store results
        angleInclination_estimates_frame(blob_index) = angleInclination_estimate;
        angleAzimuth_estimates_frame(blob_index) = angleAzimuth_estimate;
        positionX_nm_estimate_frame(blob_index) = positionX_nm_estimate;
        positionY_nm_estimate_frame(blob_index) = positionY_nm_estimate;
        defocus_estimates_frame(blob_index) = defocus_estimate;
        positionX_nm_errors_frame(blob_index) = positionX_nm_error;
        positionY_nm_errors_frame(blob_index) = positionY_nm_error;
        angleInclination_errors_frame(blob_index) = angleInclination_error;
        angleAzimuth_errors_frame(blob_index) = angleAzimuth_error;

        fprintf('    Error: (Δx, Δy) = (%i, %i)\n', round(positionX_nm_error), round(positionY_nm_error));
        fprintf('    Error: (Δθ, Δφ) = (%i°, %i°)\n', round(angleInclination_error*180/pi), round(angleAzimuth_error*180/pi));

        % draw the patch on the original image
        % if length(x_array) == 1
        %     current_color = 'red';
        % else
        %     red_value = 1 - (blob_index - 1) / (length(x_array) - 1); % Red decreases
        %     blue_value = (blob_index - 1) / (length(x_array) - 1);    % Blue increase            
        %     current_color = uint8([red_value, 0, blue_value] * 255);
        % end
        current_color = 'red';

        % Draw a red rectangle on the original image for the patch
        outlined_image = insertShape(outlined_image, 'Rectangle', ...
                                     [patch_start_x_px, patch_start_y_px, patch_width_px, patch_width_px], ...
                                     'Color', current_color, 'LineWidth', 1);
        outlined_image = insertShape(outlined_image, 'Circle', [centroid_x_px, centroid_y_px, 0], 'Color', 'red', 'LineWidth', 1);


    end % end loop over blobs

    subplot(1,3,2);
    imshow(outlined_image);
    title('Original Image with Patch Outlines');

    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    positionX_nm_estimates{frame_index} = positionX_nm_estimate_frame;
    positionY_nm_estimates{frame_index} = positionY_nm_estimate_frame;
    defocus_estimates{frame_index} = defocus_estimates_frame;
    positionX_nm_errors{frame_index} = positionX_nm_errors_frame;
    positionY_nm_errors{frame_index} = positionY_nm_errors_frame;
    angleInclination_errors{frame_index} = angleInclination_errors_frame;
    angleAzimuth_errors{frame_index} = angleAzimuth_errors_frame;


end % end loop over frames

% Flatten the array (convert cell array to a numeric array)
angleInclination_estimates = [angleInclination_estimates{:}];
angleAzimuth_estimates = [angleAzimuth_estimates{:}];
positionX_nm_estimates = [positionX_nm_estimates{:}];
positionY_nm_estimates = [positionY_nm_estimates{:}];
defocus_estimates = [defocus_estimates{:}];
positionX_nm_errors = [positionX_nm_errors{:}];
positionY_nm_errors = [positionY_nm_errors{:}];
angleInclination_errors = [angleInclination_errors{:}];
angleAzimuth_errors = [angleAzimuth_errors{:}];

% And convert to image coords for drawing later
positionX_px_estimates = [];
positionY_px_estimates = [];
for i=1:number_of_spots
    positionX_px_estimates(end+1) = nm_to_px(positionX_nm_estimates(i), pixel_size_nm, image_size_px, 'x');
    positionY_px_estimates(end+1) = nm_to_px(positionY_nm_estimates(i), pixel_size_nm, image_size_px, 'y');
end

disp('Ground truth vs estimates')
fprintf('x_tru = [%s]\n', sprintf('%.2f, ', positionX_nm_array(1:end)));
fprintf('x_est = [%s]\n', sprintf('%.2f, ', positionX_nm_estimates(1:end)));
fprintf('y_tru = [%s]\n', sprintf('%.2f, ', positionY_nm_array(1:end)));
fprintf('y_est = [%s]\n', sprintf('%.2f, ', positionY_nm_estimates(1:end)));
fprintf('inc_tru = [%s]\n', sprintf('%.2f, ', angleInclination_array(1:end)*180/pi));
fprintf('inc_est = [%s]\n', sprintf('%.2f, ', angleInclination_estimates(1:end)*180/pi));
fprintf('az_tru = [%s]\n', sprintf('%.2f, ', angleAzimuth_array(1:end)*180/pi));
fprintf('az_est = [%s]\n', sprintf('%.2f, ', angleAzimuth_estimates(1:end)*180/pi));

% Calculate errors
% mean_position_errors = (positionX_nm_errors + positionY_nm_errors) / 2;
% mean_position_error = mean(mean_position_errors);
% mean_positionX_error = mean(positionX_nm_errors);
% mean_positionY_error = mean(positionY_nm_errors);
% mean_angleInclination_error = mean(angleInclination_errors);
% mean_angleAzimuth_error = mean(angleAzimuth_errors);
% fprintf('----------\n');
% fprintf('Mean |error|: (Δx, Δy) = (%.3f, %.3f)\n', mean_positionX_error, mean_positionY_error);
% fprintf('Mean |error|: (Δθ, Δφ) = (%.3f°, %.3f°)\n', mean_angleInclination_error*180/pi, mean_angleAzimuth_error*180/pi);
% fprintf('----------\n');
% fprintf('inclination = [%s]\n', sprintf('%.2f, ', angleInclination_array(1:end)*180/pi));
% fprintf('inc_err = [%s]\n', sprintf('%.2f, ', angleInclination_errors(1:end)*180/pi));
% fprintf('loc_err = [%s]\n', sprintf('%.2f, ', mean_position_errors(1:end)));
fprintf('x_err = [%s]\n', sprintf('%.3f, ', positionX_nm_errors(1:end)));
fprintf('y_err = [%s]\n', sprintf('%.3f, ', positionY_nm_errors(1:end)));
fprintf('inc_err = [%s]\n', sprintf('%.3f, ', angleInclination_errors(1:end)*180/pi));
fprintf('az_err = [%s]\n', sprintf('%.3f, ', angleAzimuth_errors(1:end)*180/pi));


% % Format the array into a string with 4 decimal places
% formattedString_1 = sprintf('%.4f, ', angleInclination_errors*180/pi);
% formattedString_2 = sprintf('%.4f, ', positionX_nm_errors);
% formattedString_1 = ['[', formattedString_1(1:end-2), ']'];
% formattedString_2 = ['[', formattedString_2(1:end-2), ']'];
% disp(formattedString_1);
% disp(formattedString_2);

% disp(par.backgroundNoise)
% disp((mean_positionX_error + mean_positionY_error)/2)

subplot(1,3,3);
imshow(psf_total_image);
hold on;
plot(positionX_px_array, positionY_px_array, 'x', 'Color', 'r', 'MarkerSize', 5, 'LineWidth', 2);
plot(positionX_px_estimates, positionY_px_estimates, 'x', 'Color', 'g', 'MarkerSize', 5, 'LineWidth', 2);
title('localisations');
