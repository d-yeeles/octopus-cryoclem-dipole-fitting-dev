%% Fitting multiple PSFs in a single frame

% Same as looop-test/fit_multiple.m, but feeding it the known simulated
% positions rather than using thunderstorm to estimate them.

close all;
clear all;

addpath(genpath('../'));

%% ----------
%% Simulate
%% ----------

inclinations = 75*pi/180:75*pi/180;%0:pi/180:2*pi;
azimuths = 0:pi/180:2*pi;%0:0;
runs = 0:0;%0:1:25;

% Global params - these will be the same whether sim or fit

number_of_spots = 1;
scalefactor = 1;
padding = 0.15; % for avoiding edges
inner_bound = padding;
outer_bound = 1 - 2*padding;
pixel_size_nm = 51.2/scalefactor;
image_size_nm = 1000;%image_size_px*pixel_size_nm; % if arranged in NxN grid, allow 1000 nm per spot
image_size_px = roundToOdd(image_size_nm/pixel_size_nm);%roundToOdd(201);%101*scalefactor); % must be odd

% Attocube params
wavelength = 500;
objectiveFocalLength = 770;
par.nPixels = image_size_px;
par.wavelength = Length(wavelength,'nm');
par.objectiveNA = 2.17;
par.objectiveFocalLength = Length(objectiveFocalLength,'mu');
par.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
par.nDiscretizationBFP = 129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
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
    
            output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/animation_lowres/vary_az/sim_inc%i_az%i_run%i.png', round(inclination_deg), round(azimuth_deg), round(run));
    
            % Need this for checking distance from neighbours
            positionX_nm_array = [];
            positionY_nm_array = [];
            angleInclination_array = [];
            angleAzimuth_array = [];
    
            tic;
    
            for i = 1:number_of_spots

            %     min_distance_nm = 1000;
            %     valid_position = false;
            % 
            %     while ~valid_position
            % 
            %         % Generate a random position avoiding edges
            %         relative_x = inner_bound + outer_bound * rand();
            %         relative_y = inner_bound + outer_bound * rand();
            % 
            %         % Convert to nm position
            %         positionX_nm = (relative_x - 0.5) * image_size_nm;
            %         positionY_nm = (relative_y - 0.5) * image_size_nm;
            % 
            %         % Check distance from all existing spots
            %         if isempty(positionX_nm_array)
            %             valid_position = true; % First spot is always valid
            %         else
            %             distances = sqrt((positionX_nm_array - positionX_nm).^2 + ...
            %                              (positionY_nm_array - positionY_nm).^2);
            %             if all(distances >= min_distance_nm)
            %                 valid_position = true;
            %             end
            %         end
            %     end
            
                positionX_nm = 0;
                positionY_nm = 0;

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
           
        end % end loop over azimuths
    
    end % end loop over inclinations

end % end loop over runs
