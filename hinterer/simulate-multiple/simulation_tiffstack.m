%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

inclinations = 0:1*(pi/180):pi/2;%0:22.5*pi/180:pi/2;
azimuths = 0:0;%0:180*(pi/180):2*pi-180*(pi/180);
runs = 1:1;

% Global params - these will be the same whether sim or fit

number_of_spots = 9;
% scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 52;%/scalefactor;
image_size_nm = sqrt(number_of_spots)*2500;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%rous   ndToOdd(201);%101*scalefactor); % must be odd

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.nPixels = image_size_px;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 301;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
backgroundNoise = 200; % taken from looking at blank bit of example data
par.nPhotons = 4000;%1000;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images
            
output_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/10spot_testing_thunderstorm/stack_hinterer.tif';
gt_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/10spot_testing_thunderstorm/ground_truth_hinterer.csv');

frames = {};
% Ground truth
positionX_nm_array = [];
positionY_nm_array = [];
angleInclination_array = [];
angleAzimuth_array = [];
frame_number_array = [];

frame_number = 0;

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)
    
    for inclination_index = 1:length(inclinations)
    
        for azimuth_index = 1:length(azimuths)
    
            frame_number = frame_number + 1;

            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);
    
            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;
    
            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);
    
            testX_array = [];
            testY_array = [];
    
            tic;
    
            for i = 1:number_of_spots

                % Random placement
                min_distance_nm = 800;
                valid_position = false;

                while ~valid_position

                    % Generate a random position avoiding edges
                    relative_x = inner_bound + outer_bound * rand();
                    relative_y = inner_bound + outer_bound * rand();

                    % Convert to nm position
                    positionX_nm = (relative_x - 0.5) * image_size_nm;
                    positionY_nm = (relative_y - 0.5) * image_size_nm;

                    % Check distance from all existing spots
                    if isempty(testX_array)
                        valid_position = true; % First spot is always valid
                    else
                        distances = sqrt((testX_array - positionX_nm).^2 + ...
                                         (testY_array - positionY_nm).^2);
                        if all(distances >= min_distance_nm)
                            valid_position = true;
                        end
                    end

                end


                % % Regular grid
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
                %     % add some random jitter within a pixel
                %     positionX_nm = positionX_nm - pixel_size_nm/2 + rand*pixel_size_nm;
                %     positionY_nm = positionY_nm - pixel_size_nm/2 + rand*pixel_size_nm;
                % end
            


                % Just in the middle
                % positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
                % positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;

                angleInclination = inclination;%pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
                angleAzimuth = azimuth;%2*pi*rand;%0; % doesn't affect anything
            
                par.position = Length([positionX_nm positionY_nm 0], 'nm'); % This is nm away from centre of image
                par.dipole = Dipole(angleInclination, angleAzimuth);
            
                % no noise for each spot that we're adding on top of each
                % other, then put noise in the final one

                if i == number_of_spots
                    par.backgroundNoise = backgroundNoise;
                    par.shotNoise = 1;
                else 
                    par.backgroundNoise = 0;
                    % par.shotNoise = 0;
                end

                psf = PSF(par);
            
                % if first iteration, use this psf image
                % later loops just add to this image
                if i == 1
                    psf_total_image = psf.image;
                else
                    psf_total_image = psf_total_image + psf.image;
                end




                % Clip values just for display
                display_image = psf_total_image;
                display_image = double(display_image); % Convert to double for calculations
                display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
                imshow(display_image)




                % These are just used to ensure valid positions
                testX_array(end+1) = positionX_nm;
                testY_array(end+1) = positionY_nm;

                positionX_nm_array(end+1) = positionX_nm + image_size_nm/2; % modified to match coords of thunderstorm results
                positionY_nm_array(end+1) = image_size_nm/2 - positionY_nm; % modified to match coords of thunderstorm results
                angleInclination_array(end+1) = angleInclination;
                angleAzimuth_array(end+1) = angleAzimuth;
                frame_number_array(end+1) = frame_number;

            end % end loop over blobs
    
            elapsed_time = toc;
            fprintf('    Generated frame in %.2f seconds\n', elapsed_time);
    
            % % Output as png
            % psf_total_image = uint32(psf_total_image);
            % display_image = double(psf_total_image); % Convert to double for calculations
            % display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
            % imwrite(display_image, output_path);

            % Output as tif
            psf_total_image = uint32(psf_total_image);

            % Store the generated frame in the frames cell array
            frames{end+1} = psf_total_image;
    
            % Clip values just for display
            display_image = psf_total_image;
            display_image = double(display_image); % Convert to double for calculations
            display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
            imshow(display_image)

            % fprintf('Simulation output to \n %s\n', output_path);
           
        end % end loop over azimuths
    
    end % end loop over inclinations

end % end loop over runs

disp(length(frames))

% Output as tif
t = Tiff(output_path, 'w');
tagstruct.ImageLength = image_size_px;  % Set image height
tagstruct.ImageWidth = image_size_px;   % Set image width
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
tagstruct.BitsPerSample = 32;  % 32-bit per pixel
tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
tagstruct.RowsPerStrip = 16;   % Strip length for compression
tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression (optional)
tagstruct.Software = 'MATLAB';
t.setTag(tagstruct);
t.write(frames{1});  % Write the first frame to the TIFF file

% Now append the rest of the frames
for i = 2:length(frames)
    t.writeDirectory();  % Create a new directory for the next frame
    t.setTag(tagstruct);  % Set the same tags for the new frame
    t.write(frames{i});   % Write the current frame to the TIFF file
end

t.close();  % Close the file when all frames are written

% Save ground truth info
fileID = fopen(gt_output_path, 'w');
fprintf(fileID, 'frame = [%s]\n', num2str(frame_number_array, ' %d ,'));
fprintf(fileID, 'x_tru = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
fprintf(fileID, 'y_tru = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
fprintf(fileID, 'inc_tru = [%s]\n', num2str(angleInclination_array, ' %d ,'));
fprintf(fileID, 'az_tru = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
fclose(fileID);

% Combine data into a matrix
data = [frame_number_array(:), positionX_nm_array(:), positionY_nm_array(:), angleInclination_array(:), angleAzimuth_array(:)];

% Combine headers and data into a cell array
headers = {'frame', 'x_tru', 'y_tru', 'inc_tru', 'az_tru'};
data_with_headers = [headers; num2cell(data)];

% Write to CSV (this will include headers and data)
writecell(data_with_headers, gt_output_path);









%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

%% ----------
%% Simulate
%% ----------

inclinations = 0:1*(pi/180):pi/2;%0:22.5*pi/180:pi/2;
azimuths = 0:0;%0:180*(pi/180):2*pi-180*(pi/180);
runs = 1:1;

% Global params - these will be the same whether sim or fit

number_of_spots = 9;
% scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 52;%/scalefactor;
image_size_nm = sqrt(number_of_spots)*2500;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%rous   ndToOdd(201);%101*scalefactor); % must be odd

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.nPixels = image_size_px;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 301;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
par.pixelSize = Length(pixel_size_nm,'nm');
par.pixelSensitivityMask = PixelSensitivity.uniform(9);
backgroundNoise = 200; % taken from looking at blank bit of example data
par.nPhotons = 4000/2;%1000;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

output_path = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/10spot_testing_thunderstorm/stack_gaussian.tif';
gt_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/10spot_testing_thunderstorm/ground_truth_gaussian.csv');

frames = {};
% Ground truth
positionX_nm_array = [];
positionY_nm_array = [];
angleInclination_array = [];
angleAzimuth_array = [];
frame_number_array = [];

frame_number = 0;

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)

    for inclination_index = 1:length(inclinations)

        for azimuth_index = 1:length(azimuths)

            frame_number = frame_number + 1;

            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);

            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;

            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);

            testX_array = [];
            testY_array = [];

            tic;

            for i = 1:number_of_spots

                % Random placement
                min_distance_nm = 800;
                valid_position = false;

                while ~valid_position

                    % Generate a random position avoiding edges
                    relative_x = inner_bound + outer_bound * rand();
                    relative_y = inner_bound + outer_bound * rand();

                    % Convert to nm position
                    positionX_nm = (relative_x - 0.5) * image_size_nm;
                    positionY_nm = (relative_y - 0.5) * image_size_nm;

                    % Check distance from all existing spots
                    if isempty(testX_array)
                        valid_position = true; % First spot is always valid
                    else
                        distances = sqrt((testX_array - positionX_nm).^2 + ...
                                         (testY_array - positionY_nm).^2);
                        if all(distances >= min_distance_nm)
                            valid_position = true;
                        end
                    end

                end


                % % Regular grid
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
                %     % add some random jitter within a pixel
                %     positionX_nm = positionX_nm - pixel_size_nm/2 + rand*pixel_size_nm;
                %     positionY_nm = positionY_nm - pixel_size_nm/2 + rand*pixel_size_nm;
                % end



                % Just in the middle
                % positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
                % positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;

                angleInclination = inclination;%pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
                angleAzimuth = azimuth;%2*pi*rand;%0; % doesn't affect anything

                par.position = Length([positionX_nm positionY_nm 0], 'nm'); % This is nm away from centre of image
                par.dipole = Dipole(angleInclination, angleAzimuth);

                % no noise for each spot that we're adding on top of each
                % other, then put noise in the final one

                if i == number_of_spots
                    par.backgroundNoise = backgroundNoise;
                    par.shotNoise = 1;
                else 
                    par.backgroundNoise = 0;
                    % par.shotNoise = 0;
                end

                psf = PSF_gaussian(par);

                % if first iteration, use this psf image
                % later loops just add to this image
                if i == 1
                    psf_total_image = psf.image;
                else
                    psf_total_image = psf_total_image + psf.image;
                end




                % Clip values just for display
                display_image = psf_total_image;
                display_image = double(display_image); % Convert to double for calculations
                display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
                imshow(display_image)




                % These are just used to ensure valid positions
                testX_array(end+1) = positionX_nm;
                testY_array(end+1) = positionY_nm;

                positionX_nm_array(end+1) = positionX_nm + image_size_nm/2; % modified to match coords of thunderstorm results
                positionY_nm_array(end+1) = image_size_nm/2 - positionY_nm; % modified to match coords of thunderstorm results
                angleInclination_array(end+1) = angleInclination;
                angleAzimuth_array(end+1) = angleAzimuth;
                frame_number_array(end+1) = frame_number;

            end % end loop over blobs

            elapsed_time = toc;
            fprintf('    Generated frame in %.2f seconds\n', elapsed_time);

            % % Output as png
            % psf_total_image = uint32(psf_total_image);
            % display_image = double(psf_total_image); % Convert to double for calculations
            % display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
            % imwrite(display_image, output_path);

            % Output as tif
            psf_total_image = uint32(psf_total_image);

            % Store the generated frame in the frames cell array
            frames{end+1} = psf_total_image;

            % Clip values just for display
            display_image = psf_total_image;
            display_image = double(display_image); % Convert to double for calculations
            display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
            imshow(display_image)

            % fprintf('Simulation output to \n %s\n', output_path);

        end % end loop over azimuths

    end % end loop over inclinations

end % end loop over runs

disp(length(frames))

% Output as tif
t = Tiff(output_path, 'w');
tagstruct.ImageLength = image_size_px;  % Set image height
tagstruct.ImageWidth = image_size_px;   % Set image width
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
tagstruct.BitsPerSample = 32;  % 32-bit per pixel
tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
tagstruct.RowsPerStrip = 16;   % Strip length for compression
tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression (optional)
tagstruct.Software = 'MATLAB';
t.setTag(tagstruct);
t.write(frames{1});  % Write the first frame to the TIFF file

% Now append the rest of the frames
for i = 2:length(frames)
    t.writeDirectory();  % Create a new directory for the next frame
    t.setTag(tagstruct);  % Set the same tags for the new frame
    t.write(frames{i});   % Write the current frame to the TIFF file
end

t.close();  % Close the file when all frames are written

% Save ground truth info
fileID = fopen(gt_output_path, 'w');
fprintf(fileID, 'frame = [%s]\n', num2str(frame_number_array, ' %d ,'));
fprintf(fileID, 'x_tru = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
fprintf(fileID, 'y_tru = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
fprintf(fileID, 'inc_tru = [%s]\n', num2str(angleInclination_array, ' %d ,'));
fprintf(fileID, 'az_tru = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
fclose(fileID);

% Combine data into a matrix
data = [frame_number_array(:), positionX_nm_array(:), positionY_nm_array(:), angleInclination_array(:), angleAzimuth_array(:)];

% Combine headers and data into a cell array
headers = {'frame', 'x_tru', 'y_tru', 'inc_tru', 'az_tru'};
data_with_headers = [headers; num2cell(data)];

% Write to CSV (this will include headers and data)
writecell(data_with_headers, gt_output_path);

