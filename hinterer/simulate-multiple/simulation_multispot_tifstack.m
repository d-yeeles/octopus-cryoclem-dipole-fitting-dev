% %% Fitting multiple PSFs in a single frame with patch-based approach
% 
% close all;
% clear all;
% 
% addpath(genpath('../'));
% 
% %% ----------
% %% Simulate
% %% ----------
% 
% number_of_frames = 100;
% 
% % Global params - these will be the same whether sim or fit
% number_of_spots = 23;           % Number of dipoles to simulate
% scalefactor = 1;
% pixel_size_nm = 52/scalefactor;
% 
% % Global image parameters
% image_width_px = 686;
% image_height_px = 438;
% global_image_width_px = image_width_px;  % Keeping both variables for clarity
% global_image_height_px = image_height_px;  % Keeping both variables for clarity
% image_width_nm = image_width_px * pixel_size_nm;
% image_height_nm = image_height_px * pixel_size_nm;  % Fixed reference to pixel_height_nm -> pixel_size_nm
% 
% % Small patch parameters for individual dipole simulation
% patch_size_px = 19;             % Size of individual patch (must be odd)
% patch_center = floor(patch_size_px/2) + 1; % Center pixel of patch (1-indexed)
% 
% % Key fix: In MATLAB's 1-indexed system, the center pixel is:
% global_center_x = floor(global_image_width_px/2) + 1; % Center pixel (1-indexed)
% global_center_y = floor(global_image_height_px/2) + 1; % Center pixel (1-indexed)
% global_center = global_center_x; % For backward compatibility with print statements
% 
% % Padding to avoid placing dipoles too close to the edges
% edge_padding = 10;  % Minimum distance from edges (in pixels)
% 
% % Attocube params
% wavelength = 500;
% objectiveFocalLength = 770;
% par.wavelength = Length(wavelength,'nm');
% par.objectiveNA = 2.17;
% par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
% par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
% par.nDiscretizationBFP = 129;
% par.pixelSize = Length(pixel_size_nm,'nm');
% par.pixelSensitivityMask = PixelSensitivity.uniform(9);
% backgroundNoise = 0; 
% par.nPhotons = 2000;
% 
% % Create a base directory for outputs
% base_output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/sims_23spot_large';
% base_output_dir = sprintf(base_output_dir, number_of_spots);
% 
% % Create output directory if it doesn't exist
% if ~exist(base_output_dir, 'dir')
%     fprintf('Creating output directory: %s\n', base_output_dir);
%     mkdir(base_output_dir);
% end
% 
% % Initialize an array to hold all images for the stack
% stack_images = cell(number_of_frames);
% % % Initialize an array to store all metadata for each frame
% % stack_metadata = cell(length(runs), length(inclinations), length(azimuths));
% 
% % Loop over runs, inclinations, and azimuths
% for frame = 1:number_of_frames
% 
%     fprintf('Generating frame %i/%i\n', frame, number_of_frames);
% 
% 
%     % We'll still generate the individual file paths for metadata purposes
%     % single_output_path = sprintf('%s/sim_theta%03i_phi%03i_run%i.tif', base_output_dir, ceil(inclination_deg), ceil(azimuth_deg), round(run));
%     data_output_path = sprintf('%s/params_frame%06i.m', base_output_dir, round(frame));
% 
%     % Create global image canvas (initialized to zeros)
%     global_image = zeros(global_image_height_px, global_image_width_px);
% 
%     % Arrays to store position and angle data for all dipoles
%     positionX_nm_array = [];
%     positionY_nm_array = [];
%     angleInclination_array = [];
%     angleAzimuth_array = [];
% 
%     tic;
% 
%     for i = 1:number_of_spots
% 
%         inclination = (pi/2)*rand();
%         azimuth = (2*pi)*rand();
% 
%         inclination_deg = inclination*180/pi;
%         azimuth_deg = azimuth*180/pi;
% 
%         % Generate random position for this dipole (in pixel coordinates)
%         % Avoid edges by using padding
%         valid_position = false;
%         min_distance_nm = 1000; % Minimum distance between dipoles in pixels
% 
%         max_attempts = 100;
%         attempt = 0;
% 
%         while ~valid_position && attempt < max_attempts
%             attempt = attempt + 1;
% 
%             % Generate position within allowed bounds (in pixels)
%             x_px = edge_padding + randi(global_image_width_px - 2*edge_padding);
%             y_px = edge_padding + randi(global_image_height_px - 2*edge_padding);
% 
%             % Calculate where to place this patch in the global image
%             start_x = x_px - patch_center + 1;
%             start_y = y_px - patch_center + 1;
%             end_x = start_x + patch_size_px - 1;
%             end_y = start_y + patch_size_px - 1;
% 
%             % Ensure patch fits within global image bounds
%             if start_x < 1 || start_y < 1 || end_x > global_image_width_px || end_y > global_image_height_px
%                 continue; % Skip this position - patch would exceed boundaries
%             end
% 
%             % Convert to nm relative to center of global image
%             % In nm coordinates: origin is at center, +x is right, +y is up
%             % In pixel coordinates: origin is top-left, +x is right, +y is down
%             % Important: store the exact pixel position where the dipole will be placed
%             positionX_nm = (x_px - global_center_x) * pixel_size_nm;
%             positionY_nm = (global_center_y - y_px) * pixel_size_nm; % Flip y-axis
% 
%             % Check distance from all existing spots
%             if isempty(positionX_nm_array)
%                 valid_position = true; % First spot is always valid
%             else
%                 distances = sqrt((positionX_nm_array - positionX_nm).^2 + ...
%                                 (positionY_nm_array - positionY_nm).^2);
%                 if all(distances >= min_distance_nm)
%                     valid_position = true;
%                 end
%             end
%         end
% 
%         if ~valid_position
%             fprintf('Warning: Could not find valid position for dipole %d after %d attempts\n', i, max_attempts);
%             continue;
%         end
% 
%         % Store the verified position coordinates
%         positionX_nm_array(end+1) = positionX_nm;
%         positionY_nm_array(end+1) = positionY_nm;
% 
%         % Set up parameters for small patch simulation
%         % Position is set to 0,0 as we want dipole at center of small patch
%         par.nPixels = patch_size_px;
%         par.position = Length([0 0 0], 'nm'); % Centered in the small patch
% 
%         % Set dipole orientation
%         angleInclination = inclination;
%         angleAzimuth = azimuth;
% 
%         par.dipole = Dipole(angleInclination, angleAzimuth);
% 
%         % Set noise parameters for individual patch
%         if i == number_of_spots
%             par.backgroundNoise = backgroundNoise;
%             par.shotNoise = 1;
%         else 
%             par.backgroundNoise = 0;
%             par.shotNoise = 0;
%         end
% 
%         % Generate PSF for this dipole on small patch
%         psf = PSF(par);
%         patch_image = psf.image;
% 
%         % Calculate where to place this patch in the global image (AGAIN to be sure)
%         start_x = x_px - patch_center + 1;
%         start_y = y_px - patch_center + 1;
%         end_x = start_x + patch_size_px - 1;
%         end_y = start_y + patch_size_px - 1;
% 
%         % Add patch to global image
%         global_image(start_y:end_y, start_x:end_x) = global_image(start_y:end_y, start_x:end_x) + patch_image;
% 
%         % Store the dipole orientation
%         angleInclination_array(end+1) = angleInclination;
%         angleAzimuth_array(end+1) = angleAzimuth;
% 
%         % Calculate patch center for verification
%         patch_center_x = start_x + floor(patch_size_px/2);
%         patch_center_y = start_y + floor(patch_size_px/2);
% 
%         % Double-check conversion for debugging
%         test_x_px = round(positionX_nm / pixel_size_nm + global_center_x);
%         test_y_px = round(global_center_y - positionY_nm / pixel_size_nm);
% 
%         % Verify exact pixel where dipole is placed
%         if abs(test_x_px - x_px) > 1 || abs(test_y_px - y_px) > 1
%             fprintf('  WARNING: POSITION MISMATCH! Generated at (%d, %d) but would place at (%d, %d)\n', ...
%                     x_px, y_px, test_x_px, test_y_px);
%         end
%     end
% 
%     elapsed_time = toc;
%     fprintf('    %.2f seconds\n', elapsed_time);
% 
%     % Store the generated image in our stack
%     stack_images{frame} = uint32(global_image);
% 
%     % Calculate image size in nm based on the global image size
%     image_width_nm = global_image_width_px * pixel_size_nm;
%     image_height_nm = global_image_height_px * pixel_size_nm;
% 
%     % % Store metadata for this frame
%     % metadata = struct();
%     % metadata.inclination_deg = inclination_deg;
%     % metadata.azimuth_deg = azimuth_deg;
%     % metadata.run = run;
%     % metadata.number_of_spots = number_of_spots;
%     % metadata.pixel_size_nm = pixel_size_nm;
%     % metadata.image_width_nm = image_width_nm;
%     % metadata.image_height_nm = image_height_nm;
%     % metadata.image_width_px = image_width_px;
%     % metadata.image_height_px = image_height_px;
%     % metadata.patch_size_px = patch_size_px;
%     % metadata.wavelength = wavelength;
%     % metadata.objectiveNA = par.objectiveNA;
%     % metadata.objectiveFocalLength = objectiveFocalLength;
%     % metadata.refractiveIndices = par.refractiveIndices;
%     % metadata.nDiscretizationBFP = par.nDiscretizationBFP;
%     % metadata.backgroundNoise = par.backgroundNoise;
%     % metadata.nPhotons = par.nPhotons;
%     % metadata.positionX_nm_array = positionX_nm_array;
%     % metadata.positionY_nm_array = positionY_nm_array;
%     % metadata.angleInclination_array = angleInclination_array;
%     % metadata.angleAzimuth_array = angleAzimuth_array;
%     % 
%     % stack_metadata{run, inclination_index, azimuth_index} = metadata;
% 
%     % Save individual metadata file (if still needed)
%     fileID = fopen(data_output_path, 'w');
%     fprintf(fileID, '%% ground truth for frame %i\n', round(frame));
%     % fprintf(fileID, '%% settings\n');
%     % fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
%     % fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
%     % fprintf(fileID, 'image_width_nm = %.3f\n', image_width_nm);
%     % fprintf(fileID, 'image_height_nm = %.3f\n', image_height_nm);
%     % fprintf(fileID, 'image_width_px = %.3f\n', image_width_px);
%     % fprintf(fileID, 'image_height_px = %.3f\n', image_height_px);
%     % fprintf(fileID, 'patch_size_px = %i\n', patch_size_px);  % Add this to track patch size
%     % fprintf(fileID, 'wavelength = %i\n', wavelength);
%     % fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
%     % fprintf(fileID, 'objectiveFocalLength = %i\n', objectiveFocalLength);
%     % fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %g ,'));
%     % fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
%     % fprintf(fileID, 'par.backgroundNoise = %.3f\n', par.backgroundNoise);
%     % fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
%     fprintf(fileID, '%% data\n');
%     fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, '%.10f, '));
%     fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, '%.10f, '));
%     fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, '%.10f, '));
%     fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, '%.10f, '));
%     fclose(fileID);
% 
% end % end loop over frames
% 
% % At this point, all simulations are done and stored in stack_images
% % Now we need to flatten our 3D cell array and save as a TIF stack
% 
% % % Calculate total number of frames
% % total_frames = length(runs) * length(inclinations) * length(azimuths);
% % fprintf('Total frames generated: %d\n', total_frames);
% 
% % Prepare metadata for the stack
% stack_output_path = sprintf('%s/stack_all_simulations.tif', base_output_dir);
% % stack_metadata_path = sprintf('%s/stack_metadata.mat', base_output_dir);
% 
% % Save the stack as a multi-page TIF
% fprintf('Saving all frames to TIF stack: %s\n', stack_output_path);
% 
% try
%     % Create a fresh file first
%     t = Tiff(stack_output_path, 'w');
%     tagstruct.ImageLength = global_image_height_px;
%     tagstruct.ImageWidth = global_image_width_px;
%     tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.BitsPerSample = 32;
%     tagstruct.SamplesPerPixel = 1;
%     tagstruct.RowsPerStrip = 16;
%     tagstruct.Compression = Tiff.Compression.LZW;
%     tagstruct.Software = 'MATLAB';
%     tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     % Add resolution tags to avoid dimension issues
%     tagstruct.XResolution = 72;
%     tagstruct.YResolution = 72;
%     tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Inch;
% 
%     % Write the frames
%     for frame_idx = 1:number_of_frames
%         % Make sure we have a valid image with non-zero dimensions
%         if isempty(stack_images{frame_idx}) || all(size(stack_images{frame_idx}) == 0)
%             fprintf('Warning: Empty image at frame %d, skipping\n', frame_idx);
%             continue;
%         end
% 
%         if frame_idx == 1
%             % First frame write
%             t.setTag(tagstruct);
%             t.write(stack_images{frame_idx});
%         else
%             % Create a new directory for subsequent frames
%             t.writeDirectory();
%             t.setTag(tagstruct);
%             t.write(stack_images{frame_idx});
%         end
% 
%         fprintf('Added frame %d/%d to stack\n', frame_idx, number_of_frames);
%     end
% 
%     t.close();
%     fprintf('Successfully saved TIF stack with %d frames\n', number_of_frames);
% 
% catch err
%     fprintf('Error saving TIF stack: %s\n', err.message);
% end
% 
% fprintf('Processing complete.\n');
% 




%% 1 spot version below, randomised within central pixel

%% Fitting multiple PSFs in a single frame with patch-based approach
% Modified to place dipoles randomly within the central pixel

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

model = 'hinterer';
base_output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_hinterer_0';

theta_values = 0:0;%[0*pi/180, 22.5*pi/180, 45*pi/180, 67.5*pi/180, 90*pi/180];
phi_values = 0:1*pi/180:2*pi-0.00000001;
number_of_frames = length(theta_values) * length(phi_values);


% Global params - these will be the same whether sim or fit
number_of_spots = 1;           % Number of dipoles to simulate
scalefactor = 1;
pixel_size_nm = 52/scalefactor;

% Global image parameters
image_width_px = 19*scalefactor;
image_height_px = 19*scalefactor;
global_image_width_px = image_width_px;  % Keeping both variables for clarity
global_image_height_px = image_height_px;  % Keeping both variables for clarity
image_width_nm = image_width_px * pixel_size_nm;
image_height_nm = image_height_px * pixel_size_nm;

% Small patch parameters for individual dipole simulation
patch_size_px = 19*scalefactor;             % Size of individual patch (must be odd)
patch_center = floor(patch_size_px/2) + 1; % Center pixel of patch (1-indexed)

% Key fix: In MATLAB's 1-indexed system, the center pixel is:
global_center_x = floor(global_image_width_px/2) + 1; % Center pixel (1-indexed)
global_center_y = floor(global_image_height_px/2) + 1; % Center pixel (1-indexed)
global_center = global_center_x; % For backward compatibility with print statements

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

% Create output directory if it doesn't exist
if ~exist(base_output_dir, 'dir')
    fprintf('Creating output directory: %s\n', base_output_dir);
    mkdir(base_output_dir);
end

% Initialize an array to hold all images for the stack
stack_images = cell(number_of_frames, 1);  % Make sure it's a column vector

blob_counter = 0;

% Loop over frames
for theta = theta_values
    for phi = phi_values
        % Start timing for this frame
        tic;
        
        blob_counter = blob_counter + 1;
    
        fprintf('Generating frame %i/%i\n', blob_counter, length(theta_values)*length(phi_values));
    
        % We'll still generate the individual file paths for metadata purposes
        data_output_path = sprintf('%s/params_frame%06i.m', base_output_dir, round(blob_counter));
    
        % Create global image canvas (initialized to zeros)
        global_image = zeros(global_image_height_px, global_image_width_px);
    
        % Arrays to store position and angle data for all dipoles
        positionX_nm_array = [];
        positionY_nm_array = [];
        angleInclination_array = [];
        angleAzimuth_array = [];
    
        for i = 1:number_of_spots
            % % Calculate the orientation angles based on frame number
            % phi_range = 360; % 0-359 degrees
            % frames_per_complete_cycle = phi_range * num_theta;
            % cycle_frame = mod(frame - 1, frames_per_complete_cycle);
            % theta_idx = floor(cycle_frame / phi_range) + 1;
            % azimuth = mod(cycle_frame, phi_range) * pi/180;
            % inclination = theta_values(theta_idx);
    
            inclination = theta;
            azimuth = phi;

            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;
    
            % Calculate the boundaries of the central pixel in nm
            pixel_min_x_nm = -pixel_size_nm/2;
            pixel_max_x_nm = pixel_size_nm/2;
            pixel_min_y_nm = -pixel_size_nm/2;
            pixel_max_y_nm = pixel_size_nm/2;
    
            % Generate a random position within the central pixel (in nm)
            random_offset_x_nm = pixel_min_x_nm + rand() * (pixel_max_x_nm - pixel_min_x_nm);
            random_offset_y_nm = pixel_min_y_nm + rand() * (pixel_max_y_nm - pixel_min_y_nm);
    
            % The center of the global image is at (0,0) in nm coordinates
            positionX_nm = random_offset_x_nm;
            positionY_nm = random_offset_y_nm;
    
            % Calculate the exact pixel position this corresponds to
            x_px = global_center_x + positionX_nm / pixel_size_nm;
            y_px = global_center_y - positionY_nm / pixel_size_nm; % Flip y-axis
    
            % For patch placement, we need integer pixel positions
            x_px_int = round(x_px);
            y_px_int = round(y_px);
    
            % Store the position coordinates in nm
            positionX_nm_array(end+1) = positionX_nm;
            positionY_nm_array(end+1) = positionY_nm;
    
            % Calculate where to place this patch in the global image
            start_x = x_px_int - patch_center + 1;
            start_y = y_px_int - patch_center + 1;
            end_x = start_x + patch_size_px - 1;
            end_y = start_y + patch_size_px - 1;
    
            % Check if patch fits within global image bounds
            if start_x < 1 || start_y < 1 || end_x > global_image_width_px || end_y > global_image_height_px
                fprintf('Warning: Patch would exceed image boundaries. Adjusting position.\n');
                % Adjust position to ensure patch fits
                if start_x < 1
                    x_px_int = x_px_int + (1 - start_x);
                    start_x = 1;
                end
                if start_y < 1
                    y_px_int = y_px_int + (1 - start_y);
                    start_y = 1;
                end
                if end_x > global_image_width_px
                    x_px_int = x_px_int - (end_x - global_image_width_px);
                    end_x = global_image_width_px;
                end
                if end_y > global_image_height_px
                    y_px_int = y_px_int - (end_y - global_image_height_px);
                    end_y = global_image_height_px;
                end
    
                % Recalculate exact position in nm after adjustment
                positionX_nm = (x_px_int - global_center_x) * pixel_size_nm;
                positionY_nm = (global_center_y - y_px_int) * pixel_size_nm; % Flip y-axis
    
                % Update stored position coordinates
                positionX_nm_array(end) = positionX_nm;
                positionY_nm_array(end) = positionY_nm;
            end
    
            % Set up parameters for small patch simulation
            par.nPixels = patch_size_px;
            par.position = Length([0 0 0], 'nm'); % Center in patch
    
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
    
            % Add patch to global image
            global_image(start_y:end_y, start_x:end_x) = global_image(start_y:end_y, start_x:end_x) + patch_image;
    
            % Store the dipole orientation
            angleInclination_array(end+1) = angleInclination;
            angleAzimuth_array(end+1) = angleAzimuth;
    
            % For verification
            fprintf('  Dipole %d positioned at (%.2f, %.2f) nm from center\n', i, positionX_nm, positionY_nm);
            fprintf('  This corresponds to pixel position (%.2f, %.2f)\n', x_px, y_px);
        end

        % Record elapsed time for this frame - moved inside both loops
        elapsed_time = toc;
        fprintf('    %.2f seconds\n', elapsed_time);

        % Store the generated image in our stack - moved inside both loops
        stack_images{blob_counter} = uint32(global_image);

        % Calculate image size in nm based on the global image size
        image_width_nm = global_image_width_px * pixel_size_nm;
        image_height_nm = global_image_height_px * pixel_size_nm;

        % Save individual metadata file - moved inside both loops
        fileID = fopen(data_output_path, 'w');
        fprintf(fileID, '%% ground truth for frame %i\n', round(blob_counter));
        fprintf(fileID, '%% data\n');
        fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, '%.10f, '));
        fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, '%.10f, '));
        fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, '%.10f, '));
        fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, '%.10f, '));
        fclose(fileID);
    end % end phi loop
end % end theta loop

% At this point, all simulations are done and stored in stack_images
% Now we save as a TIF stack

% Prepare metadata for the stack
stack_output_path = sprintf('%s/stack_all_simulations.tif', base_output_dir);

% Save the stack as a multi-page TIF
fprintf('Saving all frames to TIF stack: %s\n', stack_output_path);

try
    % Create a fresh file first
    t = Tiff(stack_output_path, 'w');
    tagstruct.ImageLength = global_image_height_px;
    tagstruct.ImageWidth = global_image_width_px;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.Compression = Tiff.Compression.LZW;
    tagstruct.Software = 'MATLAB';
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    % Add resolution tags to avoid dimension issues
    tagstruct.XResolution = 72;
    tagstruct.YResolution = 72;
    tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Inch;

    % Write the frames
    for frame_idx = 1:number_of_frames
        % Make sure we have a valid image with non-zero dimensions
        if isempty(stack_images{frame_idx}) || all(size(stack_images{frame_idx}) == 0)
            fprintf('Warning: Empty image at frame %d, skipping\n', frame_idx);
            continue;
        end

        if frame_idx == 1
            % First frame write
            t.setTag(tagstruct);
            t.write(stack_images{frame_idx});
        else
            % Create a new directory for subsequent frames
            t.writeDirectory();
            t.setTag(tagstruct);
            t.write(stack_images{frame_idx});
        end

        fprintf('Added frame %d/%d to stack\n', frame_idx, number_of_frames);
    end

    t.close();
    fprintf('Successfully saved TIF stack with %d frames\n', number_of_frames);

catch err
    fprintf('Error saving TIF stack: %s\n', err.message);
end

fprintf('Processing complete.\n');
