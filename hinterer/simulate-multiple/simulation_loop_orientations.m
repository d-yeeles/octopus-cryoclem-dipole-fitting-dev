%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

inclinations = 0:22.5*(pi/180):pi/2;%0:22.5*pi/180:pi/2;
azimuths = 0:1*(pi/180):2*pi-1*(pi/180);%0:4*pi/180:2*pi;
runs = 1:5;

% Global params - these will be the same whether sim or fit

number_of_spots = 1;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 52/scalefactor;
image_size_nm = sqrt(number_of_spots)*1000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%rous   ndToOdd(201);%101*scalefactor); % must be odd

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
par.backgroundNoise = 200; % taken from looking at blank bit of example data
par.nPhotons = 1e9;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)
    
    for inclination_index = 1:length(inclinations)
    
        for azimuth_index = 1:length(azimuths)
    
            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);
    
            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;
    
            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);
    
            output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_finaltest/sims_highN/sim_inc%03i_az%03i_run%i.tif', round(inclination_deg), round(azimuth_deg), round(run));
            data_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_finaltest/sims_highN/params_inc%03i_az%03i_run%i.m', round(inclination_deg), round(azimuth_deg), round(run));

            % Need this for checking distance from neighbours
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];
    
            tic;
    
            for i = 1:number_of_spots
                                
                % min_distance_nm = 1000;
                % valid_position = false;
                % 
                % while ~valid_position
                % 
                %     % Generate a random position avoiding edges
                %     relative_x = inner_bound + outer_bound * rand();
                %     relative_y = inner_bound + outer_bound * rand();
                % 
                %     % Convert to nm position
                %     positionX_nm = (relative_x - 0.5) * image_size_nm;
                %     positionY_nm = (relative_y - 0.5) * image_size_nm;
                % 
                %     % Check distance from all existing spots
                %     if isempty(positionX_nm_array)
                %         valid_position = true; % First spot is always valid
                %     else
                %         distances = sqrt((positionX_nm_array - positionX_nm).^2 + ...
                %                          (positionY_nm_array - positionY_nm).^2);
                %         if all(distances >= min_distance_nm)
                %             valid_position = true;
                %         end
                %     end
                % 
                % end
            
                positionX_nm = -pixel_size_nm/2 + rand*pixel_size_nm;
                positionY_nm = -pixel_size_nm/2 + rand*pixel_size_nm;

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
    
            % % Output as png
            % psf_total_image = uint32(psf_total_image);
            % display_image = double(psf_total_image); % Convert to double for calculations
            % display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
            % imwrite(display_image, output_path);

            % Output as tif
            psf_total_image = uint32(psf_total_image);
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
    
            % Clip values just for display
            display_image = imread(output_path);
            display_image = double(display_image); % Convert to double for calculations
            display_image = (display_image - min(display_image(:))) / (max(display_image(:)) - min(display_image(:)));
            imshow(display_image)

            fprintf('Simulation output to \n %s\n', output_path);

            % Save ground truth info
    
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
