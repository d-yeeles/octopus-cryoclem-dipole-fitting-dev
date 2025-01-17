
% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

inclinations = 0:pi/180:pi/2;
azimuths = pi/4:pi/4;
runs = 0:1:25;

% Global params - these will be the same whether sim or fit

number_of_spots = 2;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = sqrt(number_of_spots)*2000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.nPixels = image_size_px;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
par.backgroundNoise = 1/number_of_spots; % taken from looking at blank bit of example data
par.nPhotons = 500;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)
    
    for inclination_index = 1:length(inclinations)
    
        for azimuth_index = 1:length(azimuths)
    
            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);
    
            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;
    
            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);
    
            output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az90/sim_inc%i_az%i_run%i.tif', round(inclination_deg), round(azimuth_deg), round(run));
    
            % Need this for checking distance from neighbours
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];
    
            tic;
    
            for i = 1:number_of_spots
                                
                min_distance_nm = 1000;
                valid_position = false;
            
                while ~valid_position
            
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
            
                angleInclination = inclination;%pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
                angleAzimuth = azimuth;%2*pi*rand;%0; % doesn't affect anything
            
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
            angleInclination_array(end+1) = angleInclination;
            angleAzimuth_array(end+1) = angleAzimuth;
    
            end % end loop over blobs
    
            elapsed_time = toc;
            fprintf('    Generated frame in %.2f seconds\n', elapsed_time);
    
            % Output as tif stack
            psf_total_image = mat2gray(psf_total_image);
            psf_total_image = uint8(255*psf_total_image);
            imwrite(psf_total_image, output_path);
            fprintf('Simulation output to \n %s\n', output_path);
    
            % Save ground truth info
            data_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az90/params_inc%i_az%i_run%i.m', round(inclination_deg), round(azimuth_deg), round(run));
    
            fileID = fopen(data_output_path, 'w');
            fprintf(fileID, '%% ground truth for sim_inc%i_az%i_run%i.tif\n', round(inclination_deg), round(azimuth_deg), round(run));
            fprintf(fileID, '%% settings\n');
            fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
            fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
            fprintf(fileID, 'image_size_nm = %.3f\n', image_size_nm);
            fprintf(fileID, 'image_size_px = %i\n', image_size_px);
            fprintf(fileID, 'wavelength = %i\n', wavelength);
            fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
            fprintf(fileID, 'objectiveFocalLength = %i\n', objectiveFocalLength);
            fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %d ,'));
            fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
            fprintf(fileID, 'par.backgroundNoise = %.3f\n', par.backgroundNoise);
            fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
            fprintf(fileID, '% data\n');
            fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
            fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
            fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, ' %d ,'));
            fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
            fclose(fileID);
           
        end % end loop over azimuths
    
    end % end loop over inclinations

end % end loop over runs















% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Simulate
%% ----------

inclinations = 0:pi/180:pi/2;
azimuths = pi/2:pi/2;
runs = 0:1:25;

% Global params - these will be the same whether sim or fit

number_of_spots = 2;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = sqrt(number_of_spots)*2000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.nPixels = image_size_px;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
par.backgroundNoise = 1/number_of_spots; % taken from looking at blank bit of example data
par.nPhotons = 500;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)
    
    for inclination_index = 1:length(inclinations)
    
        for azimuth_index = 1:length(azimuths)
    
            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);
    
            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;
    
            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);
    
            output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az180/sim_inc%i_az%i_run%i.tif', round(inclination_deg), round(azimuth_deg), round(run));
    
            % Need this for checking distance from neighbours
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];
    
            tic;
    
            for i = 1:number_of_spots
                                
                min_distance_nm = 1000;
                valid_position = false;
            
                while ~valid_position
            
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
            
                angleInclination = inclination;%pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
                angleAzimuth = azimuth;%2*pi*rand;%0; % doesn't affect anything
            
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
            angleInclination_array(end+1) = angleInclination;
            angleAzimuth_array(end+1) = angleAzimuth;
    
            end % end loop over blobs
    
            elapsed_time = toc;
            fprintf('    Generated frame in %.2f seconds\n', elapsed_time);
    
            % Output as tif stack
            psf_total_image = mat2gray(psf_total_image);
            psf_total_image = uint8(255*psf_total_image);
            imwrite(psf_total_image, output_path);
            fprintf('Simulation output to \n %s\n', output_path);
    
            % Save ground truth info
            data_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az180/params_inc%i_az%i_run%i.m', round(inclination_deg), round(azimuth_deg), round(run));
    
            fileID = fopen(data_output_path, 'w');
            fprintf(fileID, '%% ground truth for sim_inc%i_az%i_run%i.tif\n', round(inclination_deg), round(azimuth_deg), round(run));
            fprintf(fileID, '%% settings\n');
            fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
            fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
            fprintf(fileID, 'image_size_nm = %.3f\n', image_size_nm);
            fprintf(fileID, 'image_size_px = %i\n', image_size_px);
            fprintf(fileID, 'wavelength = %i\n', wavelength);
            fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
            fprintf(fileID, 'objectiveFocalLength = %i\n', objectiveFocalLength);
            fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %d ,'));
            fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
            fprintf(fileID, 'par.backgroundNoise = %.3f\n', par.backgroundNoise);
            fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
            fprintf(fileID, '% data\n');
            fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
            fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
            fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, ' %d ,'));
            fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
            fclose(fileID);
           
        end % end loop over azimuths
    
    end % end loop over inclinations

end % end loop over runs












% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Simulate
%% ----------

inclinations = 0:pi/180:pi/2;
azimuths = 3*pi/4:3*pi/4;
runs = 0:1:25;

% Global params - these will be the same whether sim or fit

number_of_spots = 2;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = sqrt(number_of_spots)*2000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.nPixels = image_size_px;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
par.backgroundNoise = 1/number_of_spots; % taken from looking at blank bit of example data
par.nPhotons = 500;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)
    
    for inclination_index = 1:length(inclinations)
    
        for azimuth_index = 1:length(azimuths)
    
            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);
    
            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;
    
            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);
    
            output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az270/sim_inc%i_az%i_run%i.tif', round(inclination_deg), round(azimuth_deg), round(run));
    
            % Need this for checking distance from neighbours
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];
    
            tic;
    
            for i = 1:number_of_spots
                                
                min_distance_nm = 1000;
                valid_position = false;
            
                while ~valid_position
            
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
            
                angleInclination = inclination;%pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
                angleAzimuth = azimuth;%2*pi*rand;%0; % doesn't affect anything
            
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
            angleInclination_array(end+1) = angleInclination;
            angleAzimuth_array(end+1) = angleAzimuth;
    
            end % end loop over blobs
    
            elapsed_time = toc;
            fprintf('    Generated frame in %.2f seconds\n', elapsed_time);
    
            % Output as tif stack
            psf_total_image = mat2gray(psf_total_image);
            psf_total_image = uint8(255*psf_total_image);
            imwrite(psf_total_image, output_path);
            fprintf('Simulation output to \n %s\n', output_path);
    
            % Save ground truth info
            data_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az270/params_inc%i_az%i_run%i.m', round(inclination_deg), round(azimuth_deg), round(run));
    
            fileID = fopen(data_output_path, 'w');
            fprintf(fileID, '%% ground truth for sim_inc%i_az%i_run%i.tif\n', round(inclination_deg), round(azimuth_deg), round(run));
            fprintf(fileID, '%% settings\n');
            fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
            fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
            fprintf(fileID, 'image_size_nm = %.3f\n', image_size_nm);
            fprintf(fileID, 'image_size_px = %i\n', image_size_px);
            fprintf(fileID, 'wavelength = %i\n', wavelength);
            fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
            fprintf(fileID, 'objectiveFocalLength = %i\n', objectiveFocalLength);
            fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %d ,'));
            fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
            fprintf(fileID, 'par.backgroundNoise = %.3f\n', par.backgroundNoise);
            fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
            fprintf(fileID, '% data\n');
            fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
            fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
            fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, ' %d ,'));
            fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
            fclose(fileID);
           
        end % end loop over azimuths
    
    end % end loop over inclinations

end % end loop over runs






















% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az90/';

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
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az90/fitting_results_hinterer.py';

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













% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az180/';

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
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az180/fitting_results_hinterer.py';

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














% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az270/';

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
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az270/fitting_results_hinterer.py';

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





























% GAUSSIAN




% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az0/';

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
        fitResult = FitPSF_gaussian(psfInit, parEst);
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az0/fitting_results_gaussian.py';

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




















% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az90/';

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
        fitResult = FitPSF_gaussian(psfInit, parEst);
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az90/fitting_results_gaussian.py';

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













% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az180/';

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
        fitResult = FitPSF_gaussian(psfInit, parEst);
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az180/fitting_results_gaussian.py';

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














% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Fit
%% ----------

% Input params
frames_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az270/';

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
        fitResult = FitPSF_gaussian(psfInit, parEst);
        elapsed_time = toc;q
        fprintf('    Time to fit: %.2f seconds\n', elapsed_time);
        angleInclination_estimate = fitResult.estimatesPositionDefocus.ML(4);
        angleAzimuth_estimate = fitResult.estimatesPositionDefocus.ML(5);

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
fitting_results_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/2spot_az270/fitting_results_gaussian.py';

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


