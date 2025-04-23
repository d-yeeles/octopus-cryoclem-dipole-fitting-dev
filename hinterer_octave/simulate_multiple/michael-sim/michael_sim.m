%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all
clear all

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

% Global params - these will be the same whether sim or fit
% (except nPixels will be the patch size when fitting)

% Attocube params
scalefactor = 1;
padding = 0.05; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;

pixel_size_nm = 50;
image_size_px = 250; % must be odd
image_size_nm = image_size_px*pixel_size_nm;
par.nPixels = image_size_px;
par.wavelength = Length(512,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(770,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
par.backgroundNoise = 2.5; % taken from looking at blank bit of example data
par.nPhotons = 4000;%1e20;

% Initialise empty position arrays for fitting later
% (PSF object wants it in nm with origin at centre, drawing stuff wants px with origin at top-left)

% !!! READ IN GROUND TRUTH LOCATIONS HERE !!!

% positionX_centre_origin_nm_array = [];
% positionY_centre_origin_nm_array = [];
positionX_topleft_origin_px_array = [];
positionY_topleft_origin_px_array = [];


data = readtable('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/michael_ground_truth5');  % Read the CSV into a table
positionX_centre_origin_nm_array = data{:, 1};  % Extract the first column into an array
positionY_centre_origin_nm_array = data{:, 2};  % Extract the second column into an array
positionX_centre_origin_nm_array = positionX_centre_origin_nm_array - image_size_nm/2;
positionY_centre_origin_nm_array = image_size_nm/2 - positionY_centre_origin_nm_array;

for i=1:length(positionX_centre_origin_nm_array)
    positionX_topleft_origin_px_array(end+1) = nm_to_px(positionX_centre_origin_nm_array(i), image_size_nm, image_size_px, 'x');
    positionY_topleft_origin_px_array(end+1) = nm_to_px(positionY_centre_origin_nm_array(i), image_size_nm, image_size_px, 'y');
end

angleInclination_array = zeros(size(positionX_centre_origin_nm_array));
angleAzimuth_array = zeros(size(positionX_centre_origin_nm_array));

number_of_spots = length(positionX_centre_origin_nm_array);


%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/';
patch_width_nm = 500;
patch_width_px = patch_width_nm/pixel_size_nm;%13; % size of NxN pixel patch around blob centroid to consider (must be odd)

% Create the centroids in image coordinates
% 'Read in' ground truth positions (this will be rpelaced by like reading in thunderstorm estimate or something)

f_array = ones(1, number_of_spots);
x_array = positionX_centre_origin_nm_array;
y_array = positionY_centre_origin_nm_array;

disp(x_array)
disp(y_array)
disp(image_size_nm)

% Loop over each frame:
files = dir(frames_dir);
valid_extensions = {'.png', '.jpg', '.jpeg', '.bmp'};%, '.tiff', '.tif'};
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
    image = imread(frame_path);

    % Convert to greyscale if not
    if size(image, 3) == 3
        image = rgb2gray(image);
    end

    % Subtract the background roughly
    image = double(image); % Convert to double for calculations
    image = image - mean(image(:)); % Subtract the mean
    image = max(min(image, 255), 0); % Clip values to [0, 255]
    image = uint8(image); % Convert back to uint8

    figure;
    subplot(2,3,1);
    image = mat2gray(image);
    image = uint8(255*image);
    % for i=1:number_of_spots
    %     psf_total_image = insertShape(psf_total_image, 'Circle', [positionX_topleft_origin_px_array(i), positionY_topleft_origin_px_array(i), 0], 'Color', 'red', 'LineWidth', 1);
    % end
    imshow(image)
    colormap(gray);
    colorbar;
    addMM=@(x) sprintf('%i',x*pixel_size_nm - image_size_nm/2);
    xticklabels(cellfun(addMM,num2cell(xticks'),'UniformOutput',false));
    yticklabels(cellfun(addMM,num2cell(yticks'),'UniformOutput',false));
    xlabel('x, nm');
    ylabel('y, nm');
    title('Simulated data');

    % Get patches around each spot and show them
    [patches, patch_centres_x_nm, patch_centres_y_nm, outlined_image] = extract_patches(f_array, x_array, y_array, frame_index, patch_width_nm, patch_width_px, image_size_nm, image_size_px, image);

    % Uncomment to make sure patches are in the right places - not
    % implemented properly, but works for a single frame
    subplot(2,3,2);
    imshow(outlined_image);
    title('Original Image with Patch Outlines');


    % Loop over each blob in a frame, run analysis

    % Initialise results arrays
    angleInclination_estimates_frame = [];
    angleAzimuth_estimates_frame = [];
    relative_positionX_centre_origin_nm_estimate_frame = [];
    relative_positionY_centre_origin_nm_estimate_frame = [];
    defocus_estimates = [];


    for blob_index = 1:length(patches)

        fprintf('Analysing blob %d/%d\n', blob_index, length(patches));



        % --------------------
        % Simulate psf
        % !!! replace this whole chunk with using the image !!!
        % --------------------
        % Not sure what FitPSF() uses to do the fitting. Does it just use the
        % psf.image part of a psf object? Or does it use the other stuff stored in that?
        % Looking at FitPSF.m, it seems to use a lot of stuff that isn't just
        % the image. Like psf.pixelSensitivityMask, psf.stageDrift, and the
        % parameter estimates. So might need to adapt the PSF simulator so that
        % it "simulates" a psf object containing the image we choose to give
        % it, which will be a patch of a frame.
        %
        % Maybe we can adapt the PSF() class/function/whatever that currently
        % generates a simulation so that it simply generates a psf object where
        % the psf.image bit = some input image. The other stuff in the psf
        % object are experimental params mostly? I think.
        %
        % If we use thunderSTORM to get an initial etimate of localisations ,we
        % can pass those as our initial estimates for this to use.
        % Oh and because the patches we consider will be centered on those
        % locations, the initial estimate will always be (0,0) by definition.
        %
        % psf.image seems to be a 2D array
        % with the values spanning 0 to like 4000 or so
        % normalise to 255 perhaps?
        % I think it's just photon counts. Currently nPhotons = 1e5

        % --- Set input parameters ---

        % Image parameters for initialised image (anything, it gets replaced)
        par.nPixels = roundToOdd(patch_width_px); % Changing this here didn't seem to affect it
        angleInclination = 0; % initial guess
        angleAzimuth = 0; % initial guess
        par.dipole = Dipole(angleInclination, angleAzimuth);
        par.position = Length([0 0 0], 'nm'); % Because initial estimate is centre of patch (check that 0,0 is centre)
        % par.nPhotons = 1e4;
        psfInit = PSF(par);

        % show only last blob as example
        if blob_index == length(patches)
            subplot(2,3,3);
            psfInit_image = mat2gray(psfInit.image);
            psfInit_image = uint8(255*psfInit_image);
            imshow(psfInit_image)
            colormap(gray);
            title('test psf before swapping');
        end

        % now swap out the psf.image for the patch image created earlier -
        % no idea if this is right or doing the right thing at all, and I
        % don't know what it should be normalised to
        patch_image = patches{blob_index};
        patch_image = double(patch_image);
        % patch_image = patch_image ./ norm(patch_image);
        % % Normalize the image to the range [0, max(psf.image(:))]
        % patch_min = min(patch_image(:));  % Find the minimum value
        % patch_max = max(patch_image(:));  % Find the maximum value
        % patch_image = (patch_image - patch_min) * (max(psf.image(:)) / (patch_max - patch_min));
        psf = psfInit;
        psf.image = patch_image; % but need to normalise it to whatever this expects...

        % show only last blob as example
        if blob_index == length(patches)
            subplot(2,3,4);
            psf_image = mat2gray(psf.image);
            psf_image = uint8(255*psf_image);
            imshow(psf_image)
            colormap(gray);
            title('test psf after swapping');
        end

        % Here, chop out psf.image and replace it with the patch image
        % created above


        % --------------------

        % Fitting
        parEst.angleInclinationEstimate = 0;%(blob_index-1)*(pi/2)/number_of_spots;%rand()*2*pi;
        parEst.angleAzimuthEstimate = 0;
        parEst.parameterBounds.x = Length([-patch_width_nm patch_width_nm], 'nm'); % Example bounds
        parEst.parameterBounds.y = Length([-patch_width_nm patch_width_nm], 'nm'); % Example bounds
        parEst.parameterBounds.defocus = Length([-10 10], 'nm'); % Example bounds
        parEst.parameterStartValues.x = Length(0, 'nm'); % Example bounds
        parEst.parameterStartValues.y = Length(0, 'nm'); % Example bounds
        parEst.parameterStartValues.defocus = Length(0, 'nm'); % Example bounds
        fitResult = FitPSF(psf, parEst);

        % Store results

        % Estimates for angles are stored in:
        %   fitResult.angleInclinationEstimate
        %   fitResult.angleAzimuthEstimate
        %
        % Position is estimated with both least squares and max likelihood
        %   fitResult.estimatesPositionDefocus.LS
        %   fitResult.estimatesPositionDefocus.ML
        %       given as [x, y, defocus]

        angleInclination_estimates_frame(blob_index) = fitResult.angleInclinationEstimate;
        angleAzimuth_estimates_frame(blob_index) = fitResult.angleAzimuthEstimate;
        relative_positionX_centre_origin_nm_estimate_frame(blob_index) = fitResult.estimatesPositionDefocus.ML(1);
        relative_positionY_centre_origin_nm_estimate_frame(blob_index) = fitResult.estimatesPositionDefocus.ML(2);
        defocus_estimates_frame(blob_index) = fitResult.estimatesPositionDefocus.ML(3);

    end % end loop over blobs

    angleInclination_estimates{frame_index} = angleInclination_estimates_frame;
    angleAzimuth_estimates{frame_index} = angleAzimuth_estimates_frame;
    relative_positionX_centre_origin_nm_estimates{frame_index} = relative_positionX_centre_origin_nm_estimate_frame;
    relative_positionY_centre_origin_nm_estimates{frame_index} = relative_positionY_centre_origin_nm_estimate_frame;
    defocus_estimates{frame_index} = defocus_estimates_frame;

end % end loop over frames

% Flatten the array (convert cell array to a numeric array)
angleInclination_estimates = [angleInclination_estimates{:}];
angleAzimuth_estimates = [angleAzimuth_estimates{:}];
relative_positionX_centre_origin_nm_estimates = [relative_positionX_centre_origin_nm_estimates{:}];
relative_positionY_centre_origin_nm_estimates = [relative_positionY_centre_origin_nm_estimates{:}];
defocus_estimates = [defocus_estimates{:}];

% Localisations were done relative to patch centres
% So to get new localisations, add new ones to those ones
positionX_centre_origin_nm_estimates = patch_centres_x_nm + relative_positionX_centre_origin_nm_estimates;
positionY_centre_origin_nm_estimates = patch_centres_y_nm + relative_positionY_centre_origin_nm_estimates;

% And convert to image coords for drawing later
positionX_topleft_origin_px_estimates = [];
positionY_topleft_origin_px_estimates = [];
for i=1:number_of_spots
    positionX_topleft_origin_px_estimates(end+1) = nm_to_px(positionX_centre_origin_nm_estimates(i), image_size_nm, image_size_px, 'x');
    positionY_topleft_origin_px_estimates(end+1) = nm_to_px(positionY_centre_origin_nm_estimates(i), image_size_nm, image_size_px, 'y');
end







% Calculate the differences
positionX_nm_error = abs(positionX_centre_origin_nm_array - positionX_centre_origin_nm_estimates);
positionY_nm_error = abs(positionY_centre_origin_nm_array - positionY_centre_origin_nm_estimates);
angleInclination_error = abs(angleInclination_array - angleInclination_estimates);
angleAzimuth_error = abs(angleAzimuth_array - angleAzimuth_estimates);
% Calculate the mean errors
mean_positionX_error = mean(positionX_nm_error);
mean_positionY_error = mean(positionY_nm_error);
mean_angleInclination_error = mean(angleInclination_error);
mean_angleAzimuth_error = mean(angleAzimuth_error);

% for i = 1:number_of_spots
%     fprintf('Ground truth: (x, y) = (%i, %i)\n', round(positionX_centre_origin_nm_array(i)), round(positionY_centre_origin_nm_array(i)));
%     fprintf('Localisation: (x, y) = (%i, %i)\n', round(positionX_centre_origin_nm_estimates(i)), round(positionY_centre_origin_nm_estimates(i)));
%     fprintf('Ground truth: (θ, φ) = (%i°, %i°)\n', round(angleInclination_array(i)*180/pi), round(angleAzimuth_array(i)*180/pi));
%     fprintf('Localisation: (θ, φ) = (%i°, %i°)\n', round(angleInclination_estimates(i)*180/pi), round(angleAzimuth_estimates(i)*180/pi));
%     fprintf('Error: (Δx, Δy) = (%i, %i)\n', round(positionX_nm_error(i)), round(positionY_nm_error(i)));
%     fprintf('Error: (Δθ, Δφ) = (%i°, %i°)\n', round(angleInclination_error(i)*180/pi), round(angleAzimuth_error(i)*180/pi));
%     fprintf('----------\n');
% end
% fprintf('Mean error: (Δx, Δy) = (%i, %i)\n', round(mean_positionX_error), round(mean_positionY_error));
% fprintf('Mean error: (Δθ, Δφ) = (%i°, %i°)\n', round(mean_angleInclination_error*180/pi), round(mean_angleAzimuth_error*180/pi));
% 
% 
% 
% errs_x = positionX_centre_origin_nm_array - positionX_centre_origin_nm_estimates;
% errs_y = positionY_centre_origin_nm_array - positionY_centre_origin_nm_estimates;
% mean_positive = (mean(errs_x(errs_x > 0)) + mean(errs_y(errs_y > 0)))/2;
% mean_negative = (mean(errs_x(errs_x < 0)) + mean(errs_y(errs_y < 0)))/2;
% fprintf('----------\n');
% fprintf('Mean error: %.2f\n', (mean(errs_x) + mean(errs_y))/2);
% fprintf('Mean abs(error): %.2f\n', (mean(positionX_nm_error) + mean(positionY_nm_error))/2);
% fprintf('Mean +ve errors: %.2f\n', mean_positive);
% fprintf('Mean -ve errors: %.2f\n', mean_negative);
% fprintf('Size of half a pixel: %.2f\n', pixel_size_nm/2);

disp('patch-relative localisations')
disp(relative_positionX_centre_origin_nm_estimates)
disp('.')
disp(relative_positionY_centre_origin_nm_estimates)
% disp('global localisations')
% disp(positionX_centre_origin_nm_estimates(1:7))
% disp('ground truth')
% disp(positionX_centre_origin_nm_array(1:7))
% disp('localisation error')
% disp(positionX_centre_origin_nm_array(1:7) - positionX_centre_origin_nm_estimates(1:7))
% disp('mean localisation error')
% disp(mean(positionX_centre_origin_nm_array(1:7) - positionX_centre_origin_nm_estimates(1:7)))
% disp('mean abs(localisation error)')
% disp(mean(abs(abs(positionX_centre_origin_nm_array(1:7)) - abs(positionX_centre_origin_nm_estimates(1:7)))))


% On original image, plot ground truth dot and localisation dot
result_image = image;
for i = 1:number_of_spots
    x_true = int32(positionX_topleft_origin_px_array(i));
    y_true = int32(positionY_topleft_origin_px_array(i));
    x_guess = int32(positionX_topleft_origin_px_estimates(i));
    y_guess = int32(positionY_topleft_origin_px_estimates(i));
    result_image = insertShape(result_image, 'Circle', [x_true, y_true, 0], 'Color', 'red', 'LineWidth', 1);
    result_image = insertShape(result_image, 'Circle', [x_guess, y_guess, 0], 'Color', 'green', 'LineWidth', 1);
end

subplot(2,3,5);
imshow(result_image);
title('localisations');
