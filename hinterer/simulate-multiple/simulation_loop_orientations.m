% %% Fitting multiple PSFs in a single frame
% 
% % Same as looop-test/fit_multiple.m, but feeding it the known simulated
% % positions rather than using thunderstorm to estimate them.
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
% model = 'mortensen';
% 
% inclinations = 67.5*(pi/180):67.5*(pi/180);%0:22.5*(pi/180):pi/2;
% azimuths = 0:12*(pi/180):2*pi-1*(pi/180);
% runs = 1:1;
% 
% % Global params - these will be the same whether sim or fit
% 
% number_of_spots = 1;
% scalefactor = 1;
% padding = 0.15; % for avoiding edges
% inner_bound = padding;
% outer_bound = 1 - 2*padding;
% pixel_size_nm = 52/scalefactor;
% image_size_nm = 988;%sqrt(number_of_spots)*1000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
% image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd
% 
% % Attocube params
% wavelength = 500;
% objectiveFocalLength = 770;
% par.nPixels = image_size_px;
% par.wavelength = Length(wavelength,'nm');
% par.objectiveNA = 2.17;
% par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
% par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
% par.nDiscretizationBFP = 129;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
% par.pixelSize = Length(pixel_size_nm,'nm');
% par.pixelSensitivityMask = PixelSensitivity.uniform(9);
% backgroundNoise = 0; % taken from looking at blank bit of example data
% par.nPhotons = 2000;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images
% 
% counter = 0;
% 
% % Loop over a bunch of random orientations, positions, photon counts
% for run = 1:length(runs)
% 
%     for inclination_index = 1:length(inclinations)
% 
%         for azimuth_index = 1:length(azimuths)
% 
%             counter = counter + 1;
% 
%             inclination = inclinations(inclination_index);
%             azimuth = azimuths(azimuth_index);
% 
%             inclination_deg = inclination*180/pi;
%             azimuth_deg = azimuth*180/pi;
% 
%             fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);
% 
%             output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_test/sim_frame%06d.tif', round(counter));
%             data_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_test/params_frame%06d.txt', round(counter));
% 
%             % Need this for checking distance from neighbours
%             positionX_nm_array = [];
%             positionY_nm_array = [];
%             angleInclination_array = [];
%             angleAzimuth_array = [];
% 
%             tic;
% 
%             for i = 1:number_of_spots
% 
%             %     % use this if want no overlaps
%             %     min_distance_nm = 1000;
%             %     valid_position = false;
%             % 
%             %     while ~valid_position
%             % 
%             %         % Generate a random position avoiding edges
%             %         relative_x = inner_bound + outer_bound * rand();
%             %         relative_y = inner_bound + outer_bound * rand();
%             % 
%             %         % Convert to nm position
%             %         positionX_nm = (relative_x - 0.5) * image_size_nm;
%             %         positionY_nm = (relative_y - 0.5) * image_size_nm;
%             % 
%             %         % Check distance from all existing spots
%             %         if isempty(positionX_nm_array)
%             %             valid_position = true; % First spot is always valid
%             %         else
%             %             distances = sqrt((positionX_nm_array - positionX_nm).^2 + ...
%             %                              (positionY_nm_array - positionY_nm).^2);
%             %             if all(distances >= min_distance_nm)
%             %                 valid_position = true;
%             %             end
%             %         end
%             % 
%             %     end
% 
%                 % use this if just want random positions
%                 positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
%                 positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
% 
%                 angleInclination = inclination;%pi*rand;%(i-1)*(pi/2)/number_of_spots;%pi*rand; % symmetric about pi/2
%                 angleAzimuth = azimuth;%2*pi*rand;%0; % doesn't affect anything
% 
%                 par.position = Length([positionX_nm positionY_nm 0], 'nm'); % This is nm away from centre of image
%                 par.dipole = Dipole(angleInclination, angleAzimuth);
% 
%                 % no noise for each spot that we're adding on top of each
%                 % other, then put noise in the final one
% 
%                 if i == number_of_spots
%                     par.backgroundNoise = backgroundNoise;
%                     par.shotNoise = 1;
%                 else 
%                     par.backgroundNoise = 0;
%                     par.shotNoise = 0;
%                 end
% 
%                 % psf = PSF(par);
%                 if strcmpi(model, 'hinterer')
%                     % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
%                     psf = PSF(par);
%                 elseif strcmpi(model, 'mortensen')
%                     % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
%                     psf = PSF_mortensen(par);
%                 else
%                     error('Unknown model type: %s', model);
%                 end
% 
%                 % if first iteration, use this psf image
%                 % later loops just add to this image
%                 if i == 1
%                     psf_total_image = psf.image;
%                 else
%                     psf_total_image = psf_total_image + psf.image;
%                 end
% 
%             positionX_nm_array(end+1) = positionX_nm;
%             positionY_nm_array(end+1) = positionY_nm;
%             angleInclination_array(end+1) = angleInclination;
%             angleAzimuth_array(end+1) = angleAzimuth;
% 
%             end % end loop over blobs
% 
%             elapsed_time = toc;
%             fprintf('    Generated frame in %.2f seconds\n', elapsed_time);
% 
%             % % Output as png
%             % psf_total_image = uint32(psf_total_image);
%             % display_image = double(psf_total_image); % Convert to double for calculations
%             % display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
%             % imwrite(display_image, output_path);
% 
%             % Output as tif
%             psf_total_image = uint32(psf_total_image);
%             t = Tiff(output_path, 'w');
%             tagstruct.ImageLength = image_size_px;  % Set image height
%             tagstruct.ImageWidth = image_size_px;   % Set image width
%             tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
%             tagstruct.BitsPerSample = 32;  % 16-bit per pixel (or 32-bit if needed)
%             tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
%             tagstruct.RowsPerStrip = 16;   % Strip length for compression
%             tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression (optional)
%             tagstruct.Software = 'MATLAB';
%             t.setTag(tagstruct);
%             t.write(psf_total_image);
%             t.close();
% 
%             % % Output as csv
%             % % writematrix(psf_total_image, output_path);
%             % writematrix(flipud(fliplr(psf_total_image)), output_path); % rotate upside-down to match Mortensen phi definition
% 
%             % % Clip values just for display
%             % display_image = imread(output_path);
%             % display_image = double(display_image); % Convert to double for calculations
%             % display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
%             % imshow(display_image)
% 
%             fprintf('Simulation output to \n %s\n', output_path);
% 
%             % Save ground truth info
% 
%             fileID = fopen(data_output_path, 'w');
%             fprintf(fileID, '%% ground truth for sim_inc%i_az%i_run%i.tif\n', round(inclination_deg), round(azimuth_deg), round(run));
%             fprintf(fileID, '%% settings\n');
%             fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
%             fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
%             fprintf(fileID, 'image_size_nm = %.3f\n', image_size_nm);
%             fprintf(fileID, 'image_size_px = %i\n', image_size_px);
%             fprintf(fileID, "par.wavelength = Length(%i,'nm')\n", wavelength);
%             fprintf(fileID, "par.pixelSize = Length(%i,'nm')\n", pixel_size_nm);
%             fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
%             fprintf(fileID, "par.objectiveFocalLength = Length(%i,'mu')\n", objectiveFocalLength);
%             fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %d ,'));
%             fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
%             fprintf(fileID, 'par.backgroundNoise = %.3f\n', par.backgroundNoise);
%             fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
%             fprintf(fileID, 'par.pixelSensitivityMask = PixelSensitivity.uniform(9)\n');
%             fprintf(fileID, '% data\n');
%             fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
%             fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
%             fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, ' %d ,'));
%             fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
%             fclose(fileID);
% 
%         end % end loop over azimuths
% 
%     end % end loop over inclinations
% 
% end % end loop over runs




%% modified to output a single ground truth / params file



%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

model = 'hinterer';

inclinations = [0*(pi/180), 22.5*(pi/180), 45*(pi/180), 67.5*(pi/180), 89*(pi/180)];
azimuths = 0:1*(pi/180):2*pi-1*(pi/180);
runs = 1:1;

% Calculate total number of frames for pre-allocation
total_frames = length(runs) * length(inclinations) * length(azimuths);

% Global params - these will be the same whether sim or fit

number_of_spots = 1;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 52/scalefactor;
image_size_nm = 988;%sqrt(number_of_spots)*1000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd

% Pre-allocate image stack
psf_stack = zeros(image_size_px, image_size_px, total_frames, 'uint32');

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
backgroundNoise = 0; % taken from looking at blank bit of example data
par.nPhotons = 2000;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

counter = 0;

% Create CSV file for all simulations
output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_hinterer_all/';
csv_output_path = [output_dir 'all_simulation_parameters.csv'];

% Create and write header to CSV file
fileID = fopen(csv_output_path, 'w');
fprintf(fileID, 'frame_number,run,inclination_deg,azimuth_deg,number_of_spots,pixel_size_nm,image_size_nm,image_size_px,wavelength,objectiveNA,objectiveFocalLength,refractiveIndices,nDiscretizationBFP,backgroundNoise,nPhotons,positionX_nm_array,positionY_nm_array,angleInclination_array,angleAzimuth_array\n');
fclose(fileID);

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)

    for inclination_index = 1:length(inclinations)

        for azimuth_index = 1:length(azimuths)

            counter = counter + 1;

            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);

            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;

            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);

            output_path = sprintf('%ssim_frame%06d.tif', output_dir, round(counter));

            % Need this for checking distance from neighbours
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];

            tic;

            for i = 1:number_of_spots

                % use this if just want random positions
                positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
                positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;

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
                    par.shotNoise = 0;
                end

                % psf = PSF(par);
                if strcmpi(model, 'hinterer')
                    % par.pixelSensitivityMask = PixelSensitivity.uniform(9);
                    psf = PSF(par);
                elseif strcmpi(model, 'mortensen')
                    % par.pixelSensitivityMask = PixelSensitivity.uniform(1);
                    psf = PSF_mortensen(par);
                else
                    error('Unknown model type: %s', model);
                end

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

            % Convert to uint32 and store in the stack
            psf_total_image = uint32(psf_total_image);
            psf_stack(:,:,counter) = psf_total_image;
            
            % Still save individual TIF if needed
            t = Tiff(output_path, 'w');
            tagstruct.ImageLength = image_size_px;  % Set image height
            tagstruct.ImageWidth = image_size_px;   % Set image width
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
            tagstruct.BitsPerSample = 32;  % 16-bit per pixel (or 32-bit if needed)
            tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
            tagstruct.RowsPerStrip = 16;   % Strip length for compression
            tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression (optional)
            tagstruct.Software = 'MATLAB';
            t.setTag(tagstruct);
            t.write(psf_total_image);
            t.close();

            fprintf('Simulation %d/%d output to \n %s\n', counter, total_frames, output_path);

            % Convert arrays to strings for CSV
            % For refractive indices
            ri_str = regexprep(mat2str(par.refractiveIndices), '\s+', ',');
            ri_str = regexprep(ri_str, '[', '');
            ri_str = regexprep(ri_str, ']', '');
            
            % For position and angle arrays, we need to handle multiple spots
            % Wrap arrays in quotes to make them valid CSV fields
            posX_str = ['"' regexprep(mat2str(positionX_nm_array), '\s+', ',') '"'];
            posY_str = ['"' regexprep(mat2str(positionY_nm_array), '\s+', ',') '"'];
            angleInc_str = ['"' regexprep(mat2str(angleInclination_array), '\s+', ',') '"'];
            angleAz_str = ['"' regexprep(mat2str(angleAzimuth_array), '\s+', ',') '"'];
            
            % Clean up the array strings (remove brackets)
            posX_str = regexprep(posX_str, '\[|\]', '');
            posY_str = regexprep(posY_str, '\[|\]', '');
            angleInc_str = regexprep(angleInc_str, '\[|\]', '');
            angleAz_str = regexprep(angleAz_str, '\[|\]', '');
            
            % Append parameters to CSV file
            fileID = fopen(csv_output_path, 'a');
            fprintf(fileID, '%d,%d,%.2f,%.2f,%d,%.3f,%.3f,%d,%d,%.2f,%d,%s,%d,%.3f,%d,%s,%s,%s,%s\n', ...
                counter, run, inclination_deg, azimuth_deg, number_of_spots, pixel_size_nm, image_size_nm, ...
                image_size_px, wavelength, par.objectiveNA, objectiveFocalLength, ri_str, ...
                par.nDiscretizationBFP, backgroundNoise, par.nPhotons, posX_str, posY_str, angleInc_str, angleAz_str);
            fclose(fileID);

        end % end loop over azimuths

    end % end loop over inclinations

end % end loop over runs

fprintf('All simulation parameters saved to:\n%s\n', csv_output_path);

%% ----------
%% Save TIFF Stack - not working right when try to open in fiji
%% ----------

% % Create the stack filename
% stack_output_path = [output_dir 'sim_stack.tif'];
% 
% % Save the stack as a multipage TIFF using Tiff class which supports uint32
% fprintf('Saving TIFF stack with %d frames...\n', total_frames);
% 
% % Open the stack file
% t = Tiff(stack_output_path, 'w');
% 
% % Define common tag structure
% tagstruct.ImageLength = image_size_px;  % Set image height
% tagstruct.ImageWidth = image_size_px;   % Set image width
% tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
% tagstruct.BitsPerSample = 32;  % 32-bit per pixel
% tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
% tagstruct.RowsPerStrip = 16;   % Strip length for compression
% tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.Software = 'MATLAB';
% 
% % Write the first image
% t.setTag(tagstruct);
% t.write(psf_stack(:,:,1));
% fprintf('  Progress: 1/%d frames written\n', total_frames);
% 
% % Write the rest of the images
% for i = 2:total_frames
%     % Create a new directory (page) for the next image
%     t.writeDirectory();
%     t.setTag(tagstruct);
%     t.write(psf_stack(:,:,i));
% 
%     % Show progress
%     if mod(i, 10) == 0 || i == total_frames
%         fprintf('  Progress: %d/%d frames written\n', i, total_frames);
%     end
% end
% 
% % Close the file
% t.close();
% 
% fprintf('TIFF stack saved to:\n%s\n', stack_output_path);