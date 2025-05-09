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




% %% Single spot images
% 
% 
% 
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
% model = 'hinterer';
% 
% inclinations = [0*(pi/180), 22.5*(pi/180), 45*(pi/180), 67.5*(pi/180), 89*(pi/180)];
% azimuths = 0:1*(pi/180):2*pi-1*(pi/180);
% runs = 1:1;
% 
% % Calculate total number of frames for pre-allocation
% total_frames = length(runs) * length(inclinations) * length(azimuths);
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
% % Pre-allocate image stack
% psf_stack = zeros(image_size_px, image_size_px, total_frames, 'uint32');
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
% % Create CSV file for all simulations
% output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_hinterer_all/';
% csv_output_path = [output_dir 'all_simulation_parameters.csv'];
% 
% % Create and write header to CSV file
% fileID = fopen(csv_output_path, 'w');
% fprintf(fileID, 'frame_number,run,inclination_deg,azimuth_deg,number_of_spots,pixel_size_nm,image_size_nm,image_size_px,wavelength,objectiveNA,objectiveFocalLength,refractiveIndices,nDiscretizationBFP,backgroundNoise,nPhotons,positionX_nm_array,positionY_nm_array,angleInclination_array,angleAzimuth_array\n');
% fclose(fileID);
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
%             output_path = sprintf('%ssim_frame%06d.tif', output_dir, round(counter));
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
%             % Convert to uint32 and store in the stack
%             psf_total_image = uint32(psf_total_image);
%             psf_stack(:,:,counter) = psf_total_image;
% 
%             % Still save individual TIF if needed
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
%             fprintf('Simulation %d/%d output to \n %s\n', counter, total_frames, output_path);
% 
%             % Convert arrays to strings for CSV
%             % For refractive indices
%             ri_str = regexprep(mat2str(par.refractiveIndices), '\s+', ',');
%             ri_str = regexprep(ri_str, '[', '');
%             ri_str = regexprep(ri_str, ']', '');
% 
%             % For position and angle arrays, we need to handle multiple spots
%             % Wrap arrays in quotes to make them valid CSV fields
%             posX_str = ['"' regexprep(mat2str(positionX_nm_array), '\s+', ',') '"'];
%             posY_str = ['"' regexprep(mat2str(positionY_nm_array), '\s+', ',') '"'];
%             angleInc_str = ['"' regexprep(mat2str(angleInclination_array), '\s+', ',') '"'];
%             angleAz_str = ['"' regexprep(mat2str(angleAzimuth_array), '\s+', ',') '"'];
% 
%             % Clean up the array strings (remove brackets)
%             posX_str = regexprep(posX_str, '\[|\]', '');
%             posY_str = regexprep(posY_str, '\[|\]', '');
%             angleInc_str = regexprep(angleInc_str, '\[|\]', '');
%             angleAz_str = regexprep(angleAz_str, '\[|\]', '');
% 
%             % Append parameters to CSV file
%             fileID = fopen(csv_output_path, 'a');
%             fprintf(fileID, '%d,%d,%.2f,%.2f,%d,%.3f,%.3f,%d,%d,%.2f,%d,%s,%d,%.3f,%d,%s,%s,%s,%s\n', ...
%                 counter, run, inclination_deg, azimuth_deg, number_of_spots, pixel_size_nm, image_size_nm, ...
%                 image_size_px, wavelength, par.objectiveNA, objectiveFocalLength, ri_str, ...
%                 par.nDiscretizationBFP, backgroundNoise, par.nPhotons, posX_str, posY_str, angleInc_str, angleAz_str);
%             fclose(fileID);
% 
%         end % end loop over azimuths
% 
%     end % end loop over inclinations
% 
% end % end loop over runs
% 
% fprintf('All simulation parameters saved to:\n%s\n', csv_output_path);
% 
% %% ----------
% %% Save TIFF Stack - not working right when try to open in fiji
% %% ----------
% 
% % % Create the stack filename
% % stack_output_path = [output_dir 'sim_stack.tif'];
% % 
% % % Save the stack as a multipage TIFF using Tiff class which supports uint32
% % fprintf('Saving TIFF stack with %d frames...\n', total_frames);
% % 
% % % Open the stack file
% % t = Tiff(stack_output_path, 'w');
% % 
% % % Define common tag structure
% % tagstruct.ImageLength = image_size_px;  % Set image height
% % tagstruct.ImageWidth = image_size_px;   % Set image width
% % tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
% % tagstruct.BitsPerSample = 32;  % 32-bit per pixel
% % tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
% % tagstruct.RowsPerStrip = 16;   % Strip length for compression
% % tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression
% % tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% % tagstruct.Software = 'MATLAB';
% % 
% % % Write the first image
% % t.setTag(tagstruct);
% % t.write(psf_stack(:,:,1));
% % fprintf('  Progress: 1/%d frames written\n', total_frames);
% % 
% % % Write the rest of the images
% % for i = 2:total_frames
% %     % Create a new directory (page) for the next image
% %     t.writeDirectory();
% %     t.setTag(tagstruct);
% %     t.write(psf_stack(:,:,i));
% % 
% %     % Show progress
% %     if mod(i, 10) == 0 || i == total_frames
% %         fprintf('  Progress: %d/%d frames written\n', i, total_frames);
% %     end
% % end
% % 
% % % Close the file
% % t.close();
% % 
% % fprintf('TIFF stack saved to:\n%s\n', stack_output_path);





%% ---------------------------------------------------------------------
%% ---------------------------------------------------------------------
%% ---------------------------------------------------------------------




% %% Multi-spot images
% 
% %% Fitting multiple PSFs with each spot simulated on its own patch then placed randomly on a larger canvas
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
% model = 'hinterer';
% 
% num_frames = 200; % Number of frames to generate
% 
% % Global params - these will be the same whether sim or fit
% number_of_spots = 9; % Total number of spots per frame
% 
% % Calculate how to distribute the spots across inclination and azimuth
% % So step size should be determined by total number of spots
% % For a square-like distribution:
% % If number_of_spots = 9, we want 3x3 grid (3 inclinations, 3 azimuths)
% % If number_of_spots = 16, we want 4x4 grid (4 inclinations, 4 azimuths)
% % grid_size = round(sqrt(number_of_spots));
% % 
% % % Create equally spaced inclination and azimuth values
% % inc_values = linspace(0, 89*(pi/180), grid_size);
% % az_values = linspace(0, 360*(pi/180), grid_size);
% % 
% % % In case number_of_spots is not a perfect square, we'll use the first number_of_spots combinations
% % total_combinations = length(inc_values) * length(az_values);
% % if total_combinations > number_of_spots
% %     fprintf('Warning: Grid size %dx%d creates %d combinations, but only %d spots requested.\n', ...
% %         length(inc_values), length(az_values), total_combinations, number_of_spots);
% %     fprintf('Only the first %d combinations will be used.\n', number_of_spots);
% % end
% 
% % Over whole stack, want to scan across all combinations
% % Total number of dipoles across all frames
% total_dipoles = number_of_spots * num_frames;
% 
% % Calculate grid size for distributing angles across the entire parameter space
% angle_grid_size = round(sqrt(total_dipoles));
% 
% % Create evenly spaced inclination and azimuth values across all frames
% all_inc_values = linspace(89.5*(pi/180), 89.5*(pi/180), angle_grid_size);%linspace(0, 89*(pi/180), angle_grid_size);
% all_az_values = linspace(0, 0, angle_grid_size);%linspace(0, 360*(pi/180), angle_grid_size);
% 
% fprintf('Created a %dx%d grid of angles (%d combinations) for %d total dipoles.\n', ...
%     length(all_inc_values), length(all_az_values), length(all_inc_values)*length(all_az_values), total_dipoles);
% 
% 
% 
% 
% 
% % Calculate total number of frames for pre-allocation
% total_frames = num_frames;
% 
% scalefactor = 1;
% pixel_size_nm = 52/scalefactor;
% 
% % Size of each individual spot patch
% patch_size_nm = 988;
% patch_size_px = roundToOdd(patch_size_nm/pixel_size_nm);
% 
% % Size of the larger canvas (square with both dimensions = 988nm*2*number_of_spots)
% canvas_size_nm = patch_size_nm * 2 * number_of_spots;
% canvas_width_nm = canvas_size_nm;
% canvas_height_nm = canvas_size_nm;
% canvas_width_px = roundToOdd(canvas_width_nm/pixel_size_nm);
% canvas_height_px = roundToOdd(canvas_height_nm/pixel_size_nm);
% 
% % Pre-allocate canvas stack
% canvas_stack = zeros(canvas_height_px, canvas_width_px, total_frames, 'uint32');
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
% % Create output directory if it doesn't exist
% output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/thunderstorm_tests/sims_hinterer_smiley/';
% if ~exist(output_dir, 'dir')
%     mkdir(output_dir);
% end
% 
% % Create CSV file for all simulations
% csv_output_path = [output_dir 'all_simulation_parameters.csv'];
% 
% % Create and write header to CSV file
% fileID = fopen(csv_output_path, 'w');
% fprintf(fileID, 'frame_index,dipole_index,inclination_deg,azimuth_deg,number_of_spots,pixel_size_nm,patch_size_nm,patch_size_px,canvas_size_nm,canvas_width_px,canvas_height_px,wavelength,objectiveNA,objectiveFocalLength,refractiveIndices,nDiscretizationBFP,backgroundNoise,nPhotons,patch_positionX_nm,patch_positionY_nm,dipole_posX_nm, dipole_posY_nm,angleInclination,angleAzimuth\n');
% fclose(fileID);
% 
% % Main loop - iterate through frames
% for frame = 1:num_frames
%     fprintf('Generating frame %d/%d\n', frame, num_frames);
% 
%     % Create empty canvas for this frame
%     canvas = zeros(canvas_height_px, canvas_width_px, 'uint32');
% 
%     % Initialize arrays to store all dipole positions for this frame (for
%     % visualisation later)
%     all_dipole_posX_nm = zeros(1, number_of_spots);
%     all_dipole_posY_nm = zeros(1, number_of_spots);
% 
%     % Output path for this frame
%     output_path = sprintf('%ssim_frame%06d.tif', output_dir, frame);
% 
%     tic;
% 
%     % % Generate dipoles for this frame
%     % dipole_counter = 0;
%     % 
%     % Create all combinations of inclination and azimuth
%     % We'll use only the first number_of_spots combinations
%     % for inc_idx = 1:length(inc_values)
%     %     for az_idx = 1:length(az_values)
%     %         dipole_counter = dipole_counter + 1;
%     % 
%     %         % Skip if we've already created enough spots
%     %         if dipole_counter > number_of_spots
%     %             break;
%     %         end
%     % 
%     %         % Get inclination and azimuth values
%     %         inclination = inc_values(inc_idx);
%     %         azimuth = az_values(az_idx);
% 
%     for dipole_idx = 1:number_of_spots
%         % Calculate the global index for this dipole
%         global_idx = (frame - 1) * number_of_spots + dipole_idx;
% 
%         % Calculate corresponding indices in the angle grid
%         [inc_idx, az_idx] = ind2sub([length(all_inc_values), length(all_az_values)], ...
%             mod(global_idx - 1, length(all_inc_values) * length(all_az_values)) + 1);
% 
%         % Get the inclination and azimuth for this dipole
%         inclination = all_inc_values(inc_idx);
%         azimuth = all_az_values(az_idx);
% 
%         inclination_deg = inclination * 180/pi;
%         azimuth_deg = azimuth * 180/pi;
% 
%         fprintf('  Processing dipole %d/%d: inc=%.2f° az=%.2f°\n', dipole_idx, number_of_spots, inclination_deg, azimuth_deg);
% 
%         % Set up parameters for this spot's patch
%         par.nPixels = patch_size_px;
% 
%         % Fixed at center of patch (0nm, 0nm)
%         positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
%         positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
% 
%         par.position = Length([positionX_nm positionY_nm 0], 'nm');
%         par.dipole = Dipole(inclination, azimuth);
% 
%         % Noise for individual patches
%         par.backgroundNoise = 0;
%         par.shotNoise = 1;
% 
%         % Generate PSF
%         if strcmpi(model, 'hinterer')
%             psf = PSF(par);
%         elseif strcmpi(model, 'mortensen')
%             psf = PSF_mortensen(par);
%         else
%             error('Unknown model type: %s', model);
%         end
% 
%         % Get the patch image
%         patch_image = psf.image;
% 
%         % Use this if you want multiple non-overlapping random placement
%         % Find a random position for this patch on the canvas
%         % Ensure patches don't overlap (minimum distance between patch centers)
%         min_distance_px = 8*patch_size_px; % Minimum distance between patch centers to avoid overlap
%         valid_position = false;
% 
%         max_attempts = 100; % Prevent infinite loops
%         attempt_count = 0;
% 
%         % Arrays to store previous patch positions for this frame
%         if dipole_idx == 1
%             % First spot - initialize arrays
%             prev_patch_pos_x_px = [];
%             prev_patch_pos_y_px = [];
%         end
% 
%         while ~valid_position && attempt_count < max_attempts
%             attempt_count = attempt_count + 1;
% 
%             % Generate a random position avoiding edges
%             max_x_px = canvas_width_px - patch_size_px;
%             max_y_px = canvas_height_px - patch_size_px;
% 
%             patch_x_px = randi([1, max_x_px]);
%             patch_y_px = randi([1, max_y_px]);
% 
%             % Calculate center of this patch
%             % patch_center_x = patch_x + patch_size_px/2;
%             % patch_center_y = patch_y + patch_size_px/2;
%             if mod(patch_size_px, 2) == 1  % Odd-sized patch
%                 % For odd patches, center is at the middle pixel
%                 patch_center_x_px = patch_x_px + (patch_size_px-1)/2;
%                 patch_center_y_px = patch_y_px + (patch_size_px-1)/2;
%             else  % Even-sized patch
%                 % For even patches, center is between pixels
%                 patch_center_x_px = patch_x_px + patch_size_px/2 - 0.5;
%                 patch_center_y_px = patch_y_px + patch_size_px/2 - 0.5;
%             end
% 
%             % Check distance from all existing patch centers
%             if isempty(prev_patch_pos_x_px)
%                 valid_position = true; % First patch is always valid
%             else
%                 % Check distances to all previous patches
%                 distances_px = sqrt((prev_patch_pos_x_px - patch_center_x_px).^2 + ...
%                                 (prev_patch_pos_y_px - patch_center_y_px).^2);
% 
%                 if all(distances_px >= min_distance_px)
%                     valid_position = true;
%                 end
%             end
%         end
% 
%         % If we couldn't find a non-overlapping position after max attempts,
%         % use the last attempted position anyway
%         if attempt_count >= max_attempts && ~valid_position
%             fprintf('    Warning: Could not find non-overlapping position for dipole %d after %d attempts\n', ...
%                    dipole_idx, max_attempts);
%         end
% 
%         % Store this patch center for future overlap checks
%         prev_patch_pos_x_px(end+1) = patch_center_x_px;
%         prev_patch_pos_y_px(end+1) = patch_center_y_px;
% 
%         % Convert to nm position (center of patch relative to canvas center)
%         % canvas_center_x_px = canvas_width_px/2;
%         % canvas_center_y_px = canvas_height_px/2;
%         % patch_posX_nm = (patch_center_x - canvas_center_x_px) * pixel_size_nm;
%         % patch_posY_nm = -(patch_center_y - canvas_center_y_px) * pixel_size_nm;
%         patch_posX_nm = px_to_nm(patch_center_x_px, pixel_size_nm, canvas_width_px, 'x');
%         patch_posY_nm = px_to_nm(patch_center_y_px, pixel_size_nm, canvas_height_px, 'y');
% 
%         dipole_posX_nm = patch_posX_nm + positionX_nm;
%         dipole_posY_nm = patch_posY_nm + positionY_nm;
% 
%         % Append to arrays to store (for visualisation later)
%         all_dipole_posX_nm(dipole_idx) = dipole_posX_nm;
%         all_dipole_posY_nm(dipole_idx) = dipole_posY_nm;
% 
%         % Place the patch on the canvas
%         % Convert patch_image to uint32 to match canvas type
%         patch_image_uint32 = uint32(patch_image);
%         canvas(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1) = ...
%             canvas(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1) + patch_image_uint32;
% 
%         % Write to CSV for each dipole
%         % Convert arrays to strings for CSV
%         refractiveIndices = ['"' regexprep(mat2str(par.refractiveIndices), '\s+', ',') '"'];
%         refractiveIndices = regexprep(refractiveIndices, '\[|\]', '');
% 
%         % Append parameters to CSV file
%         fileID = fopen(csv_output_path, 'a');
%         fprintf(fileID, '%d,%d,%.2f,%.2f,%d,%.3f,%.3f,%d,%.3f,%d,%d,%d,%.2f,%d,%s,%d,%.3f,%d,%f,%f,%f,%f,%f,%f\n', ...
%             frame, dipole_idx, inclination_deg, azimuth_deg, number_of_spots, pixel_size_nm, patch_size_nm, ...
%             patch_size_px, canvas_size_nm, canvas_width_px, canvas_height_px, wavelength, ...
%             par.objectiveNA, objectiveFocalLength, refractiveIndices, par.nDiscretizationBFP, backgroundNoise, ...
%             par.nPhotons, patch_posX_nm, patch_posY_nm, dipole_posX_nm, dipole_posY_nm, inclination, azimuth);
%         fclose(fileID);
%     end
% 
%     %     % Break outer loop too if we've created enough spots
%     %     if dipole_counter >= number_of_spots
%     %         break;
%     %     end
%     % end
% 
%     % Add background noise to the final canvas
%     if backgroundNoise > 0
%         noise = poissrnd(backgroundNoise, canvas_height_px, canvas_width_px);
%         canvas = canvas + uint32(noise);
%     end
% 
%     elapsed_time = toc;
%     fprintf('  Generated frame in %.2f seconds\n', elapsed_time);
% 
%     % Store in the stack
%     canvas_stack(:,:,frame) = canvas;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     % % -----------------------------------------------------------
%     % % VISUALISATION
%     % % -----------------------------------------------------------
%     % 
%     % % Create RGB visualization
%     % % Scale canvas to 8-bit for visualization
%     % vis_canvas = double(canvas);
%     % vis_canvas = vis_canvas / max(vis_canvas(:)) * 255;
%     % vis_canvas = uint8(vis_canvas);
%     % 
%     % % Create RGB image (grayscale base)
%     % rgb_vis = cat(3, vis_canvas, vis_canvas, vis_canvas);
%     % 
%     % % Draw red rectangles around patch boundaries
%     % for i = 1:length(prev_patch_pos_x_px)
%     %     % Calculate top-left corner of patch
%     %     patch_x_px = round(prev_patch_pos_x_px(i) - patch_size_px/2);
%     %     patch_y_px = round(prev_patch_pos_y_px(i) - patch_size_px/2);
%     % 
%     %     % Ensure coordinates are within bounds
%     %     patch_x_px = max(1, min(patch_x_px, canvas_width_px - patch_size_px));
%     %     patch_y_px = max(1, min(patch_y_px, canvas_height_px - patch_size_px));
%     % 
%     %     % Draw red rectangle (borders only)
%     %     % Top horizontal line
%     %     rgb_vis(patch_y_px, patch_x_px:patch_x_px+patch_size_px-1, 1) = 255;
%     %     rgb_vis(patch_y_px, patch_x_px:patch_x_px+patch_size_px-1, 2) = 0;
%     %     rgb_vis(patch_y_px, patch_x_px:patch_x_px+patch_size_px-1, 3) = 0;
%     % 
%     %     % Bottom horizontal line
%     %     rgb_vis(patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1, 1) = 255;
%     %     rgb_vis(patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1, 2) = 0;
%     %     rgb_vis(patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1, 3) = 0;
%     % 
%     %     % Left vertical line
%     %     rgb_vis(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px, 1) = 255;
%     %     rgb_vis(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px, 2) = 0;
%     %     rgb_vis(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px, 3) = 0;
%     % 
%     %     % Right vertical line
%     %     rgb_vis(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px+patch_size_px-1, 1) = 255;
%     %     rgb_vis(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px+patch_size_px-1, 2) = 0;
%     %     rgb_vis(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px+patch_size_px-1, 3) = 0;
%     % end
%     % 
%     % % DEBUGGING: Print current dipole position to verify values
%     % fprintf('  Dipole %d position: (%.2f, %.2f) nm\n', dipole_counter, dipole_posX_nm, dipole_posY_nm);
%     % 
%     % % DEBUGGING: Print out all collected dipole positions
%     % fprintf('  Plotting positions for %d dipoles:\n', dipole_counter);
%     % for i = 1:dipole_counter
%     %     fprintf('    Dipole %d: (%.2f, %.2f) nm\n', i, all_dipole_posX_nm(i), all_dipole_posY_nm(i));
%     % end
%     % 
%     % % Plot all dipole positions
%     % for i = 1:dipole_counter
%     %     % Skip zeros (in case there's a dipole at exactly 0,0, this would be extremely unlikely)
%     %     if all_dipole_posX_nm(i) == 0 && all_dipole_posY_nm(i) == 0 && i > 1
%     %         fprintf('    Warning: Skipping dipole %d at (0,0) nm - likely uninitialized\n', i);
%     %         continue;
%     %     end
%     % 
%     %     % Convert dipole position from nm to pixels
%     %     dipole_x_px = nm_to_px(all_dipole_posX_nm(i), pixel_size_nm, canvas_width_px, 'x');
%     %     dipole_y_px = nm_to_px(all_dipole_posY_nm(i), pixel_size_nm, canvas_height_px, 'y');
%     % 
%     %     fprintf('    Dipole %d converted to pixels: (%.2f, %.2f) px\n', i, dipole_x_px, dipole_y_px);
%     % 
%     %     % Ensure coordinates are within bounds and convert to integer
%     %     dipole_x_px = round(max(1, min(dipole_x_px, canvas_width_px)));
%     %     dipole_y_px = round(max(1, min(dipole_y_px, canvas_height_px)));
%     % 
%     %     % Draw green dot (5x5 pixel square for visibility)
%     %     dot_size = 1;
%     %     for dy = -floor(dot_size/2):floor(dot_size/2)
%     %         for dx = -floor(dot_size/2):floor(dot_size/2)
%     %             pixel_x = dipole_x_px + dx;
%     %             pixel_y = dipole_y_px + dy;
%     % 
%     %             % Check bounds
%     %             if pixel_x >= 1 && pixel_x <= canvas_width_px && ...
%     %                pixel_y >= 1 && pixel_y <= canvas_height_px
%     % 
%     %                 % Make green dot
%     %                 rgb_vis(pixel_y, pixel_x, 1) = 0;
%     %                 rgb_vis(pixel_y, pixel_x, 2) = 255;
%     %                 rgb_vis(pixel_y, pixel_x, 3) = 0;
%     %             end
%     %         end
%     %     end
%     % end
%     % 
%     % % Save the visualization as PNG
%     % vis_output_path = sprintf('%svis_frame%06d.png', output_dir, frame);
%     % imwrite(rgb_vis, vis_output_path);
%     % fprintf('  Visualization saved to %s\n', vis_output_path);
%     % 
%     % % -----------------------------------------------------------
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     % Save individual TIF
%     t = Tiff(output_path, 'w');
%     tagstruct.ImageLength = canvas_height_px;
%     tagstruct.ImageWidth = canvas_width_px;
%     tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.BitsPerSample = 32;
%     tagstruct.SamplesPerPixel = 1;
%     tagstruct.RowsPerStrip = 16;
%     tagstruct.Compression = Tiff.Compression.LZW;
%     tagstruct.Software = 'MATLAB';
%     t.setTag(tagstruct);
%     t.write(canvas);
%     t.close();
% 
%     fprintf('  Saved to %s\n', output_path);
% end % End frame loop
% 
% fprintf('All simulation parameters saved to:\n%s\n', csv_output_path);
% 
% %% Save TIFF Stack (uncomment if needed)
% % stack_output_path = [output_dir 'sim_stack.tif'];
% % t = Tiff(stack_output_path, 'w');
% % tagstruct.ImageLength = canvas_height_px;
% % tagstruct.ImageWidth = canvas_width_px;
% % tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
% % tagstruct.BitsPerSample = 32;
% % tagstruct.SamplesPerPixel = 1;
% % tagstruct.RowsPerStrip = 16;
% % tagstruct.Compression = Tiff.Compression.LZW;
% % tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% % tagstruct.Software = 'MATLAB';
% % 
% % t.setTag(tagstruct);
% % t.write(canvas_stack(:,:,1));
% % 
% % for i = 2:total_frames
% %     t.writeDirectory();
% %     t.setTag(tagstruct);
% %     t.write(canvas_stack(:,:,i));
% % end
% % 
% % t.close();
% % fprintf('TIFF stack saved to:\n%s\n', stack_output_path);




%% 40,000 frames temporary

%% Multi-spot images

%% Fitting multiple PSFs with each spot simulated on its own patch then placed randomly on a larger canvas

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

model = 'mortensen';

num_frames = 40000; % Number of frames to generate - modified to 40,000
total_dipoles_per_frame = 10; % Total number of spots per frame - modified to 10

% Define the specific angles to use
inc_values_deg = [0, 22.5, 45, 67.5, 89.5]; % Inclination values in degrees
az_values_deg = 0:22.5:337.5; % Azimuth values in degrees (0, 22.5, 45, ..., 337.5)

% Convert to radians
inc_values_rad = inc_values_deg * (pi/180);
az_values_rad = az_values_deg * (pi/180);

% Calculate total number of angle combinations
num_inc_values = length(inc_values_rad);
num_az_values = length(az_values_rad);
total_combinations = num_inc_values * num_az_values;

% Total number of dipoles across all frames
total_dipoles = total_dipoles_per_frame * num_frames;

% Ensure we can fit all combinations equally
% Calculate how many complete sets of all combinations we need
num_complete_sets = ceil(total_dipoles / total_combinations);
fprintf('Total dipoles: %d, Total angle combinations: %d\n', total_dipoles, total_combinations);
fprintf('Using %d complete sets of all combinations\n', num_complete_sets);

% Create arrays to store all angle combinations
all_incs = zeros(1, total_dipoles);
all_azs = zeros(1, total_dipoles);

% Fill the arrays with repeated sets of all combinations
dipole_idx = 1;
for set_idx = 1:num_complete_sets
    % If we've generated enough dipoles, break
    if dipole_idx > total_dipoles
        break;
    end
    
    % Generate one complete set of all combinations
    for inc_idx = 1:num_inc_values
        for az_idx = 1:num_az_values
            % If we've generated enough dipoles, break
            if dipole_idx > total_dipoles
                break;
            end
            
            all_incs(dipole_idx) = inc_values_rad(inc_idx);
            all_azs(dipole_idx) = az_values_rad(az_idx);
            dipole_idx = dipole_idx + 1;
        end
        if dipole_idx > total_dipoles
            break;
        end
    end
end

% Shuffle the arrays to randomize the order
rng('shuffle'); % Set random seed based on current time
shuffle_idx = randperm(length(all_incs));
all_incs = all_incs(shuffle_idx);
all_azs = all_azs(shuffle_idx);

% Parameters for simulation
scalefactor = 1;
pixel_size_nm = 52/scalefactor;

% Size of each individual spot patch
patch_size_nm = 988;
patch_size_px = roundToOdd(patch_size_nm/pixel_size_nm);

% Size of the larger canvas
canvas_size_nm = patch_size_nm * 2 * total_dipoles_per_frame;
canvas_width_nm = canvas_size_nm;
canvas_height_nm = canvas_size_nm;
canvas_width_px = roundToOdd(canvas_width_nm/pixel_size_nm);
canvas_height_px = roundToOdd(canvas_height_nm/pixel_size_nm);

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
output_dir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/thunderstorm_tests/sims_mortensen_40k/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create CSV file for all simulations
csv_output_path = [output_dir 'all_simulation_parameters_1_5000.csv'];

% Create and write header to CSV file
fileID = fopen(csv_output_path, 'w');
fprintf(fileID, 'frame_index,dipole_index,inclination_deg,azimuth_deg,number_of_spots,pixel_size_nm,patch_size_nm,patch_size_px,canvas_size_nm,canvas_width_px,canvas_height_px,wavelength,objectiveNA,objectiveFocalLength,refractiveIndices,nDiscretizationBFP,backgroundNoise,nPhotons,patch_positionX_nm,patch_positionY_nm,dipole_posX_nm,dipole_posY_nm,angleInclination,angleAzimuth\n');
fclose(fileID);

% Main loop - iterate through frames
for frame = 1:5000
    if mod(frame, 100) == 0 || frame == 1 || frame == num_frames
        fprintf('Generating frame %d/%d\n', frame, num_frames);
    end
    
    % Create empty canvas for this frame
    canvas = zeros(canvas_height_px, canvas_width_px, 'uint32');

    % Initialize arrays to store all dipole positions for this frame
    all_dipole_posX_nm = zeros(1, total_dipoles_per_frame);
    all_dipole_posY_nm = zeros(1, total_dipoles_per_frame);

    % Output path for this frame
    output_path = sprintf('%ssim_frame%06d.tif', output_dir, frame);
    
    tic;
    
    % Arrays to store previous patch positions for this frame
    prev_patch_pos_x_px = [];
    prev_patch_pos_y_px = [];
    
    % Generate dipoles for this frame
    for dipole_idx = 1:total_dipoles_per_frame
        % Calculate the global index for this dipole
        global_idx = (frame - 1) * total_dipoles_per_frame + dipole_idx;
        
        % Get the inclination and azimuth for this dipole
        inclination = all_incs(global_idx);
        azimuth = all_azs(global_idx);
            
        inclination_deg = inclination * 180/pi;
        azimuth_deg = azimuth * 180/pi;
        
        if mod(frame, 100) == 0 || frame == 1 || frame == num_frames
            fprintf('  Processing dipole %d/%d: inc=%.2f° az=%.2f°\n', dipole_idx, total_dipoles_per_frame, inclination_deg, azimuth_deg);
        end
        
        % Set up parameters for this spot's patch
        par.nPixels = patch_size_px;
        
        % Fixed at center of patch with small random offset
        positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
        positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
        
        par.position = Length([positionX_nm positionY_nm 0], 'nm');
        par.dipole = Dipole(inclination, azimuth);
        
        % Noise for individual patches
        par.backgroundNoise = 0;
        par.shotNoise = 1;
        
        % Generate PSF
        if strcmpi(model, 'hinterer')
            psf = PSF(par);
        elseif strcmpi(model, 'mortensen')
            psf = PSF_mortensen(par);
        else
            error('Unknown model type: %s', model);
        end
        
        % Get the patch image
        patch_image = psf.image;
        
        % Find a random position for this patch on the canvas
        min_distance_px = 4*patch_size_px; % Minimum distance between patch centers
        valid_position = false;

        max_attempts = 10000; 
        attempt_count = 0;
        
        % Calculate the edge avoidance margin in pixels (2000nm)
        edge_margin_nm = 2000;
        edge_margin_px = round(edge_margin_nm / pixel_size_nm);
        
        while ~valid_position && attempt_count < max_attempts
            attempt_count = attempt_count + 1;

            % Generate a random position avoiding edges by 2000nm
            min_x_px = edge_margin_px + 1;
            min_y_px = edge_margin_px + 1;
            max_x_px = canvas_width_px - patch_size_px - edge_margin_px;
            max_y_px = canvas_height_px - patch_size_px - edge_margin_px;
            
            % Ensure we have a valid region (in case canvas is too small)
            if max_x_px <= min_x_px || max_y_px <= min_y_px
                error('Canvas too small to accommodate patches with edge margin. Consider increasing canvas size or reducing edge margin.');
            end

            patch_x_px = randi([min_x_px, max_x_px]);
            patch_y_px = randi([min_y_px, max_y_px]);

            % Calculate center of this patch
            if mod(patch_size_px, 2) == 1  % Odd-sized patch
                patch_center_x_px = patch_x_px + (patch_size_px-1)/2;
                patch_center_y_px = patch_y_px + (patch_size_px-1)/2;
            else  % Even-sized patch
                patch_center_x_px = patch_x_px + patch_size_px/2 - 0.5;
                patch_center_y_px = patch_y_px + patch_size_px/2 - 0.5;
            end

            % Check distance from all existing patch centers
            if isempty(prev_patch_pos_x_px)
                valid_position = true; % First patch is always valid
            else
                % Check distances to all previous patches
                distances_px = sqrt((prev_patch_pos_x_px - patch_center_x_px).^2 + ...
                                (prev_patch_pos_y_px - patch_center_y_px).^2);

                if all(distances_px >= min_distance_px)
                    valid_position = true;
                end
            end
        end

        % If we couldn't find a non-overlapping position, use the last attempted position
        if attempt_count >= max_attempts && ~valid_position && (mod(frame, 100) == 0 || frame == 1 || frame == num_frames)
            fprintf('    Warning: Could not find non-overlapping position for dipole %d after %d attempts\n', ...
                   dipole_idx, max_attempts);
        end
        
        % Store this patch center for future overlap checks
        prev_patch_pos_x_px(end+1) = patch_center_x_px;
        prev_patch_pos_y_px(end+1) = patch_center_y_px;
        
        % Convert to nm position
        patch_posX_nm = px_to_nm(patch_center_x_px, pixel_size_nm, canvas_width_px, 'x');
        patch_posY_nm = px_to_nm(patch_center_y_px, pixel_size_nm, canvas_height_px, 'y');
        
        dipole_posX_nm = patch_posX_nm + positionX_nm;
        dipole_posY_nm = patch_posY_nm + positionY_nm;

        % Store dipole positions
        all_dipole_posX_nm(dipole_idx) = dipole_posX_nm;
        all_dipole_posY_nm(dipole_idx) = dipole_posY_nm;

        % Place the patch on the canvas
        patch_image_uint32 = uint32(patch_image);
        canvas(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1) = ...
            canvas(patch_y_px:patch_y_px+patch_size_px-1, patch_x_px:patch_x_px+patch_size_px-1) + patch_image_uint32;
        
        % Write to CSV for each dipole
        % Convert arrays to strings for CSV
        refractiveIndices = ['"' regexprep(mat2str(par.refractiveIndices), '\s+', ',') '"'];
        refractiveIndices = regexprep(refractiveIndices, '\[|\]', '');
        
        % % Append parameters to CSV file (only for a sample of frames to avoid huge CSV)
        % if frame <= 100 || mod(frame, 1000) == 0 || frame == num_frames
        fileID = fopen(csv_output_path, 'a');
        fprintf(fileID, '%d,%d,%.2f,%.2f,%d,%.3f,%.3f,%d,%.3f,%d,%d,%d,%.2f,%d,%s,%d,%.3f,%d,%f,%f,%f,%f,%f,%f\n', ...
            frame, dipole_idx, inclination_deg, azimuth_deg, total_dipoles_per_frame, pixel_size_nm, patch_size_nm, ...
            patch_size_px, canvas_size_nm, canvas_width_px, canvas_height_px, wavelength, ...
            par.objectiveNA, objectiveFocalLength, refractiveIndices, par.nDiscretizationBFP, backgroundNoise, ...
            par.nPhotons, patch_posX_nm, patch_posY_nm, dipole_posX_nm, dipole_posY_nm, inclination, azimuth);
        fclose(fileID);
        % end
    end
    
    % Add background noise to the final canvas if needed
    if backgroundNoise > 0
        noise = poissrnd(backgroundNoise, canvas_height_px, canvas_width_px);
        canvas = canvas + uint32(noise);
    end
    
    elapsed_time = toc;
    if mod(frame, 100) == 0 || frame == 1 || frame == num_frames
        fprintf('  Generated frame in %.2f seconds\n', elapsed_time);
    end
    
    % Save individual TIF
    t = Tiff(output_path, 'w');
    tagstruct.ImageLength = canvas_height_px;
    tagstruct.ImageWidth = canvas_width_px;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.Compression = Tiff.Compression.LZW;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);
    t.write(canvas);
    t.close();
    
    if mod(frame, 100) == 0 || frame == 1 || frame == num_frames
        fprintf('  Saved to %s\n', output_path);
    end
    
    % Print progress every 1000 frames
    if mod(frame, 1000) == 0 && frame > 0
        fprintf('Progress: %.2f%% (%d/%d frames completed)\n', frame/num_frames*100, frame, num_frames);
    end
end

fprintf('All simulation parameters saved to:\n%s\n', csv_output_path);

fprintf('Simulation complete: %d frames with %d dipoles per frame\n', num_frames, total_dipoles_per_frame);
fprintf('Used %d angle combinations (%d inclinations × %d azimuths)\n', total_combinations, num_inc_values, num_az_values);