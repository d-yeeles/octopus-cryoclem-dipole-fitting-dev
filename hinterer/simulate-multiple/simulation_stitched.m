%% Fitting multiple PSFs in a 3x3 grid

% Modified version that creates 9 separate dipole images and arranges them in a 3x3 grid
% The output params file will contain the absolute positions within the combined image

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

inclinations = 0:22.5*(pi/180):pi/2;
azimuths = 0:45*(pi/180):2*pi-1*(pi/180);
runs = 1:1;

% Global params - these will be the same whether sim or fit

number_of_spots = 9;
grid_size = 3; % 3x3 grid
scalefactor = 1;
padding = 0.15; % for avoiding edges

% Parameters for each individual dipole image
single_image_size_px = roundToOdd(19); % Size of each individual dipole image (must be odd)
pixel_size_nm = 52/scalefactor;
single_image_size_nm = single_image_size_px * pixel_size_nm;

% Size of the final combined image
combined_image_size_px = single_image_size_px * grid_size;
combined_image_size_nm = single_image_size_nm * grid_size;

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
backgroundNoise = 0; % taken from looking at blank bit of example data
par.nPhotons = 2000; % number of photons per spot

% Loop over a bunch of random orientations, positions, photon counts
for run = 1:length(runs)
    for inclination_index = 1:length(inclinations)
        for azimuth_index = 1:length(azimuths)

            inclination = inclinations(inclination_index);
            azimuth = azimuths(azimuth_index);

            inclination_deg = inclination*180/pi;
            azimuth_deg = azimuth*180/pi;

            fprintf('Running inc=%.2f az=%.2f\n', inclination_deg, azimuth_deg);

            output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/sims_9spot/stitched/sim_theta%03i_phi%03i_run%i.tif', ceil(inclination_deg), ceil(azimuth_deg), round(run));
            data_output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/patching_tests/sims_9spot/stitched/params_theta%03i_phi%03i_run%i.m', ceil(inclination_deg), ceil(azimuth_deg), round(run));

            % Arrays to store global positions and angles
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];

            % Create a combined image to hold all dipole images
            combined_image = zeros(combined_image_size_px, combined_image_size_px);

            tic;

            % Generate a 3x3 grid of dipole images
            for row = 1:grid_size
                for col = 1:grid_size
                    spot_index = (row-1)*grid_size + col;
                    
                    % Calculate the top-left corner of this grid cell in the combined image
                    grid_start_x = (col-1) * single_image_size_px + 1;
                    grid_start_y = (row-1) * single_image_size_px + 1;
                    
                    % Calculate the center of this grid cell in nm, relative to the combined image center
                    center_x_nm = ((col-1) + 0.5) * single_image_size_nm - combined_image_size_nm/2;
                    center_y_nm = ((row-1) + 0.5) * single_image_size_nm - combined_image_size_nm/2;
                    
                    % Add small random offset within the pixel (just like in original)
                    random_offset_x = -pixel_size_nm/2 + rand*pixel_size_nm;
                    random_offset_y = -pixel_size_nm/2 + rand*pixel_size_nm;
                    
                    % Final position in nm, relative to the combined image center
                    positionX_nm = center_x_nm + random_offset_x;
                    positionY_nm = center_y_nm + random_offset_y;
                    
                    % Use the same inclination and azimuth for all dipoles (as in original)
                    angleInclination = inclination;
                    angleAzimuth = azimuth;
                    
                    % Set parameters for this dipole image
                    par.nPixels = single_image_size_px;
                    par.position = Length([0 0 0], 'nm'); % Center the dipole in its own image
                    par.dipole = Dipole(angleInclination, angleAzimuth);
                    
                    % Set noise parameters - we'll add noise to the final combined image
                    par.backgroundNoise = 0;
                    par.shotNoise = 0;
                    
                    % Generate the dipole image
                    psf = PSF(par);
                    
                    % Add this dipole image to the combined image at the appropriate position
                    combined_image(grid_start_y:grid_start_y+single_image_size_px-1, ...
                                  grid_start_x:grid_start_x+single_image_size_px-1) = ...
                                  combined_image(grid_start_y:grid_start_y+single_image_size_px-1, ...
                                                grid_start_x:grid_start_x+single_image_size_px-1) + psf.image;
                    
                    % Store the global position and angles
                    positionX_nm_array(end+1) = positionX_nm;
                    positionY_nm_array(end+1) = positionY_nm;
                    angleInclination_array(end+1) = angleInclination;
                    angleAzimuth_array(end+1) = angleAzimuth;
                end
            end
            
            % Now add noise to the combined image
            % Add background noise
            if backgroundNoise > 0
                noise = backgroundNoise * randn(size(combined_image));
                combined_image = combined_image + noise;
            end
            
            % Add shot noise (Poisson) if enabled
            % Note: Assuming the original intention was to add shot noise to the final image
            combined_image = poissrnd(combined_image);
            
            % Ensure non-negative values
            combined_image = max(0, combined_image);

            elapsed_time = toc;
            fprintf('    Generated frame in %.2f seconds\n', elapsed_time);

            % Output as tif
            combined_image = uint32(combined_image);
            t = Tiff(output_path, 'w');
            tagstruct.ImageLength = combined_image_size_px;  % Set image height
            tagstruct.ImageWidth = combined_image_size_px;   % Set image width
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
            tagstruct.BitsPerSample = 32;  % 32-bit per pixel
            tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
            tagstruct.RowsPerStrip = 16;   % Strip length for compression
            tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression (optional)
            tagstruct.Software = 'MATLAB';
            t.setTag(tagstruct);
            t.write(combined_image);
            t.close();

            fprintf('Simulation output to \n %s\n', output_path);

            % Save ground truth info with global positions
            fileID = fopen(data_output_path, 'w');
            fprintf(fileID, '%% ground truth for sim_inc%i_az%i_run%i.tif\n', round(inclination_deg), round(azimuth_deg), round(run));
            fprintf(fileID, '%% settings\n');
            fprintf(fileID, 'number_of_spots = %i\n', number_of_spots);
            fprintf(fileID, 'grid_size = %i\n', grid_size);
            fprintf(fileID, 'pixel_size_nm = %.3f\n', pixel_size_nm);
            fprintf(fileID, 'image_size_px = %i\n', combined_image_size_px);
            fprintf(fileID, 'image_size_nm = %.3f\n', combined_image_size_nm);
            fprintf(fileID, 'wavelength = %i\n', wavelength);
            fprintf(fileID, 'par.objectiveNA = %.2f\n', par.objectiveNA);
            fprintf(fileID, 'objectiveFocalLength = %i\n', objectiveFocalLength);
            fprintf(fileID, 'par.refractiveIndices = [%s]\n', num2str(par.refractiveIndices, ' %d ,'));
            fprintf(fileID, 'par.nDiscretizationBFP = %i\n', par.nDiscretizationBFP);
            fprintf(fileID, 'par.backgroundNoise = %.3f\n', backgroundNoise);
            fprintf(fileID, 'par.nPhotons = %i\n', par.nPhotons);
            fprintf(fileID, '%% data - positions are relative to center of combined image\n');
            fprintf(fileID, 'positionX_nm_array = [%s]\n', num2str(positionX_nm_array, ' %d ,'));
            fprintf(fileID, 'positionY_nm_array = [%s]\n', num2str(positionY_nm_array, ' %d ,'));
            fprintf(fileID, 'angleInclination_array = [%s]\n', num2str(angleInclination_array, ' %d ,'));
            fprintf(fileID, 'angleAzimuth_array = [%s]\n', num2str(angleAzimuth_array, ' %d ,'));
            fclose(fileID);

        end % end loop over azimuths
    end % end loop over inclinations
end % end loop over runs