%% Fitting multiple PSFs in a single frame with patch-based approach
% Modified to output a single TIF stack and sequentially numbered frames

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

model = 'mortensen';

inclinations = 45*pi/180:45*pi/180;%45*pi/180:22.5*(pi/180):pi/2 + 1*(pi/180); % Inclination angles
azimuths = 0:0;%0:45*(pi/180):2*pi - 1*(pi/180);    % Azimuth angles
runs = 1:100;

% Global params - these will be the same whether sim or fit
number_of_spots = 1;           % Number of dipoles to simulate
scalefactor = 1;
pixel_size_nm = 52/scalefactor;

% Global image parameters
image_width_px = 19;
image_height_px = 19;
global_image_width_px = image_width_px;  % Keeping both variables for clarity
global_image_height_px = image_height_px;  % Keeping both variables for clarity
image_width_nm = image_width_px * pixel_size_nm;
image_height_nm = image_height_px * pixel_size_nm;

% Small patch parameters for individual dipole simulation
patch_size_px = 19;             % Size of individual patch (must be odd)
patch_center = floor(patch_size_px/2) + 1; % Center pixel of patch (1-indexed)

% Key fix: In MATLAB's 1-indexed system, the center pixel is:
global_center_x = floor(global_image_width_px/2) + 1; % Center pixel (1-indexed)
global_center_y = floor(global_image_height_px/2) + 1; % Center pixel (1-indexed)

% Padding to avoid placing dipoles too close to the edges
edge_padding = 0;  % Minimum distance from edges (in pixels)

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129;
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
backgroundNoise = 0; 
par.nPhotons = 2000;

% Setup output directory
output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_45';
if ~exist(output_dir, 'dir')
    fprintf('Creating output directory: %s\n', output_dir);
    mkdir(output_dir);
end

% Prepare stack output path
stack_output_path = fullfile(output_dir, 'sim_stack.tif');

% Frame counter for sequential numbering
frame_counter = 1;
total_frames = length(runs) * length(inclinations) * length(azimuths);

% Array to store all images for the stack
all_frames = cell(total_frames, 1);

% Loop over runs, inclinations, and azimuths
for run = 1:length(runs)
    for inclination_index = 1:length(inclinations)
        for azimuth_index = 1:length(azimuths)
            % Get current inclination and azimuth for reference
            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);
            
            % For randomized angles, uncomment these:
            % inclination = 2*pi*rand();
            % azimuth = (pi/2)*rand();

            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;

            % Create sequential numbered output path
            individual_output_path = fullfile(output_dir, sprintf('sim_%06d.tif', frame_counter));
            data_output_path = fullfile(output_dir, sprintf('params_%06d.m', frame_counter));

            fprintf('Generating frame %d of %d (inc=%.2f az=%.2f)\n', frame_counter, total_frames, inclination_deg, azimuth_deg);

            % Create global image canvas (initialized to zeros)
            global_image = zeros(global_image_height_px, global_image_width_px);

            % Arrays to store position and angle data for all dipoles
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];

            tic;

            % Debug: Display center point
            fprintf('Global image center pixel: (%d, %d)\n', global_center_x, global_center_y);

            if number_of_spots == 1
                % Fix for single spot case
                i = 1;  % Loop counter for the single spot
                
                % For each spot, use the same inclination with a small adjustment in azimuth
                inclination = inclinations(inclination_index);
                azimuth = azimuths(azimuth_index);

                inclination_deg = inclination*180/pi;
                azimuth_deg = azimuth*180/pi;

                fprintf('  Spot %d: inc=%.2f az=%.2f\n', i, inclination_deg, azimuth_deg);

                % Position the dipole in the center of the image with a small random offset
                x_px = global_center_x;  % Center of the global image
                y_px = global_center_y;  % Center of the global image

                positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
                positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;

                positionX_nm_array(end+1) = positionX_nm;
                positionY_nm_array(end+1) = positionY_nm;

                % Set up parameters for small patch simulation
                % Position is set to 0,0 as we want dipole at center of small patch
                par.nPixels = patch_size_px;
                par.position = Length([0 0 0], 'nm'); % Centered in the small patch

                % Set dipole orientation
                angleInclination = inclination;
                angleAzimuth = azimuth;

                par.dipole = Dipole(angleInclination, angleAzimuth);

                % Set noise parameters for individual patch
                if i == number_of_spots
                    par.backgroundNoise = backgroundNoise;
                    par.shotNoise = 1;
                else 
                    par.backgroundNoise = 0;
                    par.shotNoise = 0;
                end

                % Generate PSF for this dipole on small patch
                if strcmpi(model, 'hinterer')
                    psf = PSF(par);
                elseif strcmpi(model, 'mortensen')
                    psf = PSF_mortensen(par);
                else
                    error('Unknown model type: %s', model);
                end

                patch_image = psf.image;

                % Calculate where to place this patch in the global image
                start_x = x_px - patch_center + 1;
                start_y = y_px - patch_center + 1;
                end_x = start_x + patch_size_px - 1;
                end_y = start_y + patch_size_px - 1;

                % Add patch to global image
                global_image(start_y:end_y, start_x:end_x) = global_image(start_y:end_y, start_x:end_x) + patch_image;

                % Store the dipole orientation
                angleInclination_array(end+1) = angleInclination;
                angleAzimuth_array(end+1) = angleAzimuth;

                % Calculate patch center for verification
                patch_center_x = start_x + floor(patch_size_px/2);
                patch_center_y = start_y + floor(patch_size_px/2);

                fprintf('    Added dipole %d at position (%.2f, %.2f) nm\n', i, positionX_nm, positionY_nm);
                fprintf('    Dipole %d patch center pixel: (%d, %d)\n', i, patch_center_x, patch_center_y);
            else
                for i = 1:number_of_spots
                    % For each spot, use the same inclination with a small adjustment in azimuth
                    inclination = inclinations(inclination_index);
                    azimuth = azimuths(azimuth_index) + i*pi/180; % each spot will be one degree apart
    
                    inclination_deg = inclination*180/pi;
                    azimuth_deg = azimuth*180/pi;
    
                    fprintf('  Spot %d: inc=%.2f az=%.2f\n', i, inclination_deg, azimuth_deg);
    
                    % Generate random position for this dipole (in pixel coordinates)
                    % Avoid edges by using padding
                    valid_position = false;
                    min_distance_px = 20; % Minimum distance between dipoles in pixels
    
                    max_attempts = 100;
                    attempt = 0;
    
                    while ~valid_position && attempt < max_attempts
                        attempt = attempt + 1;
    
                        % Generate position within allowed bounds (in pixels)
                        x_px = edge_padding + randi(global_image_width_px - 2*edge_padding);
                        y_px = edge_padding + randi(global_image_height_px - 2*edge_padding);
    
                        % Calculate where to place this patch in the global image
                        start_x = x_px - patch_center + 1;
                        start_y = y_px - patch_center + 1;
                        end_x = start_x + patch_size_px - 1;
                        end_y = start_y + patch_size_px - 1;
    
                        % Ensure patch fits within global image bounds
                        if start_x < 1 || start_y < 1 || end_x > global_image_width_px || end_y > global_image_height_px
                            continue; % Skip this position - patch would exceed boundaries
                        end
    
                        % Convert to nm relative to center of global image
                        % In nm coordinates: origin is at center, +x is right, +y is up
                        % In pixel coordinates: origin is top-left, +x is right, +y is down
                        % Important: store the exact pixel position where the dipole will be placed
                        positionX_nm = (x_px - global_center_x) * pixel_size_nm;
                        positionY_nm = (global_center_y - y_px) * pixel_size_nm; % Flip y-axis
    
                        % Check distance from all existing spots
                        if isempty(positionX_nm_array)
                            valid_position = true; % First spot is always valid
                        else
                            distances = sqrt((positionX_nm_array - positionX_nm).^2 + ...
                                            (positionY_nm_array - positionY_nm).^2);
                            if all(distances >= min_distance_px * pixel_size_nm)
                                valid_position = true;
                            end
                        end
                    end
    
                    if ~valid_position
                        fprintf('Warning: Could not find valid position for dipole %d after %d attempts\n', i, max_attempts);
                        continue;
                    end
    
                    % Store the verified position coordinates
                    positionX_nm_array(end+1) = positionX_nm;
                    positionY_nm_array(end+1) = positionY_nm;
    
                    % Set up parameters for small patch simulation
                    % Position is set to 0,0 as we want dipole at center of small patch
                    par.nPixels = patch_size_px;
                    par.position = Length([0 0 0], 'nm'); % Centered in the small patch
    
                    % Set dipole orientation
                    angleInclination = inclination;
                    angleAzimuth = azimuth;
    
                    par.dipole = Dipole(angleInclination, angleAzimuth);
    
                    % Set noise parameters for individual patch
                    if i == number_of_spots
                        par.backgroundNoise = backgroundNoise;
                        par.shotNoise = 1;
                    else 
                        par.backgroundNoise = 0;
                        par.shotNoise = 0;
                    end
    
                    % Generate PSF for this dipole on small patch
                    if strcmpi(model, 'hinterer')
                        psf = PSF(par);
                    elseif strcmpi(model, 'mortensen')
                        psf = PSF_mortensen(par);
                    else
                        error('Unknown model type: %s', model);
                    end
    
                    patch_image = psf.image;
    
                    % Calculate where to place this patch in the global image
                    start_x = x_px - patch_center + 1;
                    start_y = y_px - patch_center + 1;
                    end_x = start_x + patch_size_px - 1;
                    end_y = start_y + patch_size_px - 1;
    
                    % Add patch to global image
                    global_image(start_y:end_y, start_x:end_x) = global_image(start_y:end_y, start_x:end_x) + patch_image;
    
                    % Store the dipole orientation
                    angleInclination_array(end+1) = angleInclination;
                    angleAzimuth_array(end+1) = angleAzimuth;
    
                    % Calculate patch center for verification
                    patch_center_x = start_x + floor(patch_size_px/2);
                    patch_center_y = start_y + floor(patch_size_px/2);
    
                    fprintf('    Added dipole %d at position (%.2f, %.2f) nm\n', i, positionX_nm, positionY_nm);
                    fprintf('    Dipole %d patch center pixel: (%d, %d)\n', i, patch_center_x, patch_center_y);
                end
            end

            elapsed_time = toc;
            fprintf('    Generated frame %d with %d dipoles in %.2f seconds\n', ...
                    frame_counter, length(positionX_nm_array), elapsed_time);
            
            % Save the individual frame using ImageJ-compatible approach
            try
                % Store the original double precision version for the stack
                all_frames{frame_counter} = global_image;
                
                % For individual files, save as 16-bit TIFF (ImageJ compatible)
                % Normalize and convert to uint16
                image_to_save = double(global_image);
                if max(image_to_save(:)) > 0
                    image_to_save = image_to_save / max(image_to_save(:)) * 65535;
                end
                image_to_save = uint16(image_to_save);
                
                % Save using simple imwrite (more ImageJ-friendly)
                imwrite(image_to_save, individual_output_path, 'tif', 'Compression', 'none');
                fprintf('    Saved individual frame to %s\n', individual_output_path);
            catch ME
                fprintf('Error saving TIFF file: %s\n', ME.message);
                
                % Try the traditional TIFF method as backup
                try
                    % Fallback to using Tiff class with no compression
                    t = Tiff(individual_output_path, 'w');
                    tagstruct.ImageLength = global_image_height_px;
                    tagstruct.ImageWidth = global_image_width_px;
                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                    tagstruct.BitsPerSample = 16;
                    tagstruct.SamplesPerPixel = 1;
                    tagstruct.Compression = Tiff.Compression.None;
                    tagstruct.Software = 'MATLAB';
                    t.setTag(tagstruct);
                    
                    % Normalize and convert to uint16
                    image_to_save = double(global_image);
                    if max(image_to_save(:)) > 0
                        image_to_save = image_to_save / max(image_to_save(:)) * 65535;
                    end
                    t.write(uint16(image_to_save));
                    t.close();
                    fprintf('    Saved individual frame using fallback method\n');
                catch ME2
                    fprintf('Error with fallback TIFF method: %s\n', ME2.message);
                    
                    % Final attempt: save as PNG
                    png_output_path = strrep(individual_output_path, '.tif', '.png');
                    try
                        % Normalize for display
                        display_image = double(global_image);
                        display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
                        imwrite(display_image, png_output_path);
                        fprintf('Saved image as PNG instead: %s\n', png_output_path);
                    catch ME3
                        fprintf('Failed to save image in any format\n');
                    end
                end
            end

            % Calculate image size in nm based on the global image size
            image_width_nm = global_image_width_px * pixel_size_nm;
            image_height_nm = global_image_height_px * pixel_size_nm;

            % Save ground truth info with frame number
            fileID = fopen(data_output_path, 'w');
            fprintf(fileID, '%% ground truth for frame %06d\n', frame_counter);
            fprintf(fileID, '%% Original parameters: inc=%i°, az=%i°, run=%i\n', ...
                   round(inclination_deg), round(azimuth_deg), run);
            fprintf(fileID, '%% settings\n');
            fprintf(fileID, 'frame_number = %i\n', frame_counter);
            fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
            fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
            fprintf(fileID, 'image_width_nm = %.3f\n', image_width_nm);
            fprintf(fileID, 'image_height_nm = %.3f\n', image_height_nm);
            fprintf(fileID, 'image_width_px = %.3f\n', image_width_px);
            fprintf(fileID, 'image_height_px = %.3f\n', image_height_px);
            fprintf(fileID, 'patch_size_px = %i\n', patch_size_px);
            fprintf(fileID, 'wavelength = %i\n', wavelength);
            fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
            fprintf(fileID, 'objectiveFocalLength = %i\n', objectiveFocalLength);
            fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %d ,'));
            fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
            fprintf(fileID, 'par.backgroundNoise = %.3f\n', par.backgroundNoise);
            fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
            fprintf(fileID, '%% data\n');
            fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, '%.10f, '));
            fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, '%.10f, '));
            fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, '%.10f, '));
            fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, '%.10f, '));
            fclose(fileID);
            
            % Increment frame counter
            frame_counter = frame_counter + 1;
        end % end loop over azimuths
    end % end loop over inclinations
end % end loop over runs

%% Create the stack file
fprintf('Creating stack file with %d frames...\n', total_frames);

try
    % Method 1: Use ImageJ-friendly approach with Bio-Formats style
    % Convert all frames to 16-bit for better compatibility
    for i = 1:total_frames
        % Normalize and convert to uint16 (most compatible with ImageJ)
        frame = double(all_frames{i});
        if max(frame(:)) > 0
            frame = frame / max(frame(:)) * 65535;
        end
        all_frames{i} = uint16(frame);
    end
    
    % Use simpler TIFF writing without compression for the stack
    % This is more compatible with ImageJ
    first_frame = all_frames{1};
    imwrite(first_frame, stack_output_path, 'tif', 'Compression', 'none');
    
    % Append additional frames
    for i = 2:total_frames
        imwrite(all_frames{i}, stack_output_path, 'tif', 'Compression', 'none', 'WriteMode', 'append');
    end
    
    fprintf('Successfully created stack file: %s\n', stack_output_path);
catch ME
    fprintf('Error creating stack file: %s\n', ME.message);
    
    % Try alternative method - use lower-level TIFF functionality
    try
        fprintf('Attempting alternative method for creating stack...\n');
        
        % Open the stack TIFF file
        t = Tiff(stack_output_path, 'w');
        
        % Use simplified tag structure with no compression
        tagstruct = struct('ImageLength', global_image_height_px, ...
                          'ImageWidth', global_image_width_px, ...
                          'Photometric', Tiff.Photometric.MinIsBlack, ...
                          'BitsPerSample', 16, ...
                          'SamplesPerPixel', 1, ...
                          'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky, ...
                          'Compression', Tiff.Compression.None, ...
                          'Software', 'MATLAB');
        
        % Convert frame to uint16 if needed
        first_frame = double(all_frames{1});
        if max(first_frame(:)) > 0
            first_frame = first_frame / max(first_frame(:)) * 65535;
        end
        first_frame = uint16(first_frame);
        
        % Write the first frame
        t.setTag(tagstruct);
        t.write(first_frame);
        
        % Write remaining frames with no compression
        for i = 2:total_frames
            t.writeDirectory();
            t.setTag(tagstruct);
            
            % Convert to uint16
            frame = double(all_frames{i});
            if max(frame(:)) > 0
                frame = frame / max(frame(:)) * 65535;
            end
            t.write(uint16(frame));
        end
        
        t.close();
        fprintf('Alternative stack creation completed: %s\n', stack_output_path);
    catch ME2
        fprintf('Both stack creation methods failed.\n');
        fprintf('Error 1: %s\n', ME.message);
        fprintf('Error 2: %s\n', ME2.message);
        
        % Final fallback: Create separate files only
        fprintf('Fallback: Individual frames were saved successfully.\n');
    end
end

fprintf('Simulation complete! Generated %d individual frames and stack.\n', total_frames);