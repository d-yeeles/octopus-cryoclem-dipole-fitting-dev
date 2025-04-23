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

number_of_simulations = 3000;
number_of_spots = 1;
scalefactor = 1;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = 1000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
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
par.backgroundNoise = 0; % taken from looking at blank bit of example data
par.nPhotons = 30*scalefactor^2;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Initialise empty position arrays for fitting later
% (PSF object wants it in nm with origin at centre, drawing stuff wants px with origin at top-left)
positionX_nm_array = [];
positionY_nm_array = [];
positionX_px_array = [];
positionY_px_array = [];
angleInclination_array = [];
angleAzimuth_array = [];

% Loop over a bunch of random orientations, positions, photon counts
for i = 1:number_of_simulations

    % Convert to nm position
    positionX_nm = -pixel_size_nm/2 + pixel_size_nm * rand;
    positionY_nm = -pixel_size_nm/2 + pixel_size_nm * rand;

    % we'll need its pixel equivalent for drawing later on
    positionX_px = nm_to_px(positionX_nm, pixel_size_nm, image_size_px, 'x');
    positionY_px = nm_to_px(positionY_nm, pixel_size_nm, image_size_px, 'y');

    angleInclination = rand * (pi/2);%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
    angleAzimuth = 0;%rand * (2*pi);%0; % doesn't affect anything

    par.position = Length([positionX_nm positionY_nm 0], 'nm'); % This is nm away from centre of image
    par.dipole = Dipole(angleInclination, angleAzimuth);

    psf = PSF(par);

    % if first iteration, use this psf image
    % later loops just add to this image
    psf_total_image = psf.image;

    positionX_nm_array(end+1) = positionX_nm;
    positionY_nm_array(end+1) = positionY_nm;
    positionX_px_array(end+1) = positionX_px;
    positionY_px_array(end+1) = positionY_px;
    angleInclination_array(end+1) = angleInclination;
    angleAzimuth_array(end+1) = angleAzimuth;
    
    % Output as tif stack
    inclination_deg = angleInclination*180/pi;
    azimuth_deg = angleAzimuth*180/pi;
    output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/testing_bias/low_res/sim_inc%i_az%i_run%i.png', round(inclination_deg), round(azimuth_deg), round(i));
    imwrite(psf_total_image, output_path);
    % fprintf('Simulation output to \n %s\n', output_path);
    
end



%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/testing_bias/low_res';
patch_width_nm = image_size_nm/2; % size of patch around blob to consider
patch_width_px = roundToOdd(patch_width_nm/pixel_size_nm);% must be odd

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

    % fprintf('----------\n');
    % fprintf('FRAME %d\n', frame_index);
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

    for blob_index = 1:number_of_spots

        % fprintf(' ∘ Blob %d/%d\n', blob_index, number_of_spots);

        % Generate throwaway PSF object, later replace the image with our masked image
        par.position = Length([0 0 0], 'nm');
        psfInit = PSF(par);
        psfInit.image = psf_image;

        % Give patch as localisation bounds
        parEst.parameterBounds.x = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.inclination = [0 pi/2];
        parEst.parameterBounds.azimuth = [0 2*pi];
        % Give centroid as initial position estimate
        parEst.parameterStartValues.x = Length(0, 'nm');
        parEst.parameterStartValues.y = Length(0, 'nm');
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
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
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

        % fprintf('    Error: (Δx, Δy) = (%i, %i)\n', round(positionX_nm_error), round(positionY_nm_error));
        % fprintf('    Error: (Δθ, Δφ) = (%i°, %i°)\n', round(angleInclination_error*180/pi), round(angleAzimuth_error*180/pi));

    end % end loop over blobs

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

% fprintf('x_tru = [%s]\n', sprintf('%.3f, ', positionX_nm_array(1:end)));
% fprintf('y_tru = [%s]\n', sprintf('%.3f, ', positionY_nm_array(1:end)));
% fprintf('x_est = [%s]\n', sprintf('%.3f, ', positionX_nm_estimates(1:end)));
% fprintf('y_est = [%s]\n', sprintf('%.3f, ', positionY_nm_estimates(1:end)));
% fprintf('x_err = [%s]\n', sprintf('%.3f, ', positionX_nm_errors(1:end)));
% fprintf('y_err = [%s]\n', sprintf('%.3f, ', positionY_nm_errors(1:end)));

fprintf('mean(x_err) = %s\n', sprintf('%.3f, ', mean(positionX_nm_errors)));
fprintf('mean(y_err) = %s\n', sprintf('%.3f, ', mean(positionY_nm_errors)));

scatter(positionX_nm_array, positionY_nm_array, 'r');
hold on;
scatter(positionX_nm_estimates, positionY_nm_estimates, 'g');
xlabel('x');
ylabel('y');
xlim([-pixel_size_nm/2 pixel_size_nm/2]);
ylim([-pixel_size_nm/2 pixel_size_nm/2]);
legend('True', 'Est');

saveas(gcf, '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/testing_bias/low_res_results.png');

disp('done low res')





close all;
clear all;











%% ----------
%% Simulate
%% ----------

% Global params - these will be the same whether sim or fit

number_of_simulations = 1500;
number_of_spots = 1;
scalefactor = 8;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = 1000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
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
par.backgroundNoise = 0; % taken from looking at blank bit of example data
par.nPhotons = 30*scalefactor^2;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Initialise empty position arrays for fitting later
% (PSF object wants it in nm with origin at centre, drawing stuff wants px with origin at top-left)
positionX_nm_array = [];
positionY_nm_array = [];
positionX_px_array = [];
positionY_px_array = [];
angleInclination_array = [];
angleAzimuth_array = [];

% Loop over a bunch of random orientations, positions, photon counts
for i = 1:number_of_simulations

    % Convert to nm position
    positionX_nm = -pixel_size_nm/2 + pixel_size_nm * rand;
    positionY_nm = -pixel_size_nm/2 + pixel_size_nm * rand;

    % we'll need its pixel equivalent for drawing later on
    positionX_px = nm_to_px(positionX_nm, pixel_size_nm, image_size_px, 'x');
    positionY_px = nm_to_px(positionY_nm, pixel_size_nm, image_size_px, 'y');

    angleInclination = rand * (pi/2);%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
    angleAzimuth = 0;%rand * (2*pi);%0; % doesn't affect anything

    par.position = Length([positionX_nm positionY_nm 0], 'nm'); % This is nm away from centre of image
    par.dipole = Dipole(angleInclination, angleAzimuth);

    psf = PSF(par);

    % if first iteration, use this psf image
    % later loops just add to this image
    psf_total_image = psf.image;

    positionX_nm_array(end+1) = positionX_nm;
    positionY_nm_array(end+1) = positionY_nm;
    positionX_px_array(end+1) = positionX_px;
    positionY_px_array(end+1) = positionY_px;
    angleInclination_array(end+1) = angleInclination;
    angleAzimuth_array(end+1) = angleAzimuth;
    
    % Output as tif stack
    inclination_deg = angleInclination*180/pi;
    azimuth_deg = angleAzimuth*180/pi;
    output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/testing_bias/high_res/sim_inc%i_az%i_run%i.png', round(inclination_deg), round(azimuth_deg), round(i));
    imwrite(psf_total_image, output_path);
    % fprintf('Simulation output to \n %s\n', output_path);
    
end



%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/testing_bias/high_res';
patch_width_nm = image_size_nm/2; % size of patch around blob to consider
patch_width_px = roundToOdd(patch_width_nm/pixel_size_nm);% must be odd

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

    % fprintf('----------\n');
    % fprintf('FRAME %d\n', frame_index);
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

    for blob_index = 1:number_of_spots

        % fprintf(' ∘ Blob %d/%d\n', blob_index, number_of_spots);

        % Generate throwaway PSF object, later replace the image with our masked image
        par.position = Length([0 0 0], 'nm');
        psfInit = PSF(par);
        psfInit.image = psf_image;

        % Give patch as localisation bounds
        parEst.parameterBounds.x = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.y = Length([-image_size_nm/2 image_size_nm/2], 'nm');
        parEst.parameterBounds.defocus = Length([-10 10], 'nm');
        parEst.parameterBounds.inclination = [0 pi/2];
        parEst.parameterBounds.azimuth = [0 2*pi];
        % Give centroid as initial position estimate
        parEst.parameterStartValues.x = Length(0, 'nm');
        parEst.parameterStartValues.y = Length(0, 'nm');
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
        % fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
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

        % fprintf('    Error: (Δx, Δy) = (%i, %i)\n', round(positionX_nm_error), round(positionY_nm_error));
        % fprintf('    Error: (Δθ, Δφ) = (%i°, %i°)\n', round(angleInclination_error*180/pi), round(angleAzimuth_error*180/pi));

    end % end loop over blobs

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

% fprintf('x_tru = [%s]\n', sprintf('%.3f, ', positionX_nm_array(1:end)));
% fprintf('y_tru = [%s]\n', sprintf('%.3f, ', positionY_nm_array(1:end)));
% fprintf('x_est = [%s]\n', sprintf('%.3f, ', positionX_nm_estimates(1:end)));
% fprintf('y_est = [%s]\n', sprintf('%.3f, ', positionY_nm_estimates(1:end)));
% fprintf('x_err = [%s]\n', sprintf('%.3f, ', positionX_nm_errors(1:end)));
% fprintf('y_err = [%s]\n', sprintf('%.3f, ', positionY_nm_errors(1:end)));

fprintf('mean(x_err) = %s\n', sprintf('%.3f, ', mean(positionX_nm_errors)));
fprintf('mean(y_err) = %s\n', sprintf('%.3f, ', mean(positionY_nm_errors)));

scatter(positionX_nm_array, positionY_nm_array, 'r');
hold on;
scatter(positionX_nm_estimates, positionY_nm_estimates, 'g');
xlabel('x');
ylabel('y');
xlim([-pixel_size_nm/2 pixel_size_nm/2]);
ylim([-pixel_size_nm/2 pixel_size_nm/2]);
legend('True', 'Est');

saveas(gcf, '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/testing_bias/high_res_results.png');

disp('done high res')