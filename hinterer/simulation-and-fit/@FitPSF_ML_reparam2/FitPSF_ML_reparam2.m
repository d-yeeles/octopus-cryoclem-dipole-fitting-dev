classdef FitPSF_ML_reparam2
    
    properties
        psf % PSF  % dave apr 2025 - removing this to avoid type coercion
        image
        % angleInclinationEstimate (1,1) {mustBeInFullRadialRange} % dave jan 2025 - commented out for adding angle optimiser
        % angleAzimuthEstimate (1,1) {mustBeInFullRadialRange} % dave jan 2025 - commented out for adding angle optimiser
        noiseEstimate (1,1) {mustBeNonnegative, mustBeGreaterThanOrEqual(noiseEstimate, 1e-5)} = 1e-5  % >= 1e-5 for numeric stability of log(psf)
        nPhotonEstimate (1,1) {mustBeNonnegative} 
        stageDrift StageDrift = NoStageDrift()

        pixelSensitivityMask = PixelSensitivity.uniform(3)
        
        parameterBounds = struct( ...
            'x', Length([-800 800], 'nm'), ...
            'y', Length([-800 800], 'nm'), ...
            'defocus', Length([-2000 2000], 'nm'), ...
            'newangle1', [-1, 1], ...       % dave jan 2025 - adding angle optimiser
            'newangle2', [-1, 1], ...       % dave jan 2025 - adding angle optimiser
            'newangle3', [0, 1], ...          % dave jan 2025 - adding angle optimiser
            'photons', [1, 1e10]);          % dave jan 2025 - adding photon count optimiser
        
        parameterStartValues = struct( ...
            'x', Length(-100 + 200 * rand(), 'nm'), ...
            'y', Length(-100 + 200 * rand(), 'nm'), ...
            'defocus', Length(-500 + 1000 * rand(), 'nm'), ...
            'newangle1', 2*(rand()-0.5), ...   % dave jan 2025 - adding angle optimiser
            'newangle2', 2*(rand()-0.5), ...  % dave jan 2025 - adding angle optimiser
            'newangle3', rand(), ...  % dave jan 2025 - adding angle optimiser
            'photons', 1000);          % dave jan 2025 - adding photon count optimiser
        
        % Fit result
        estimatesPositionDefocus
        
        % BFP type indicator
        model = 'hinterer'  % Default to Hinterer BackFocalPlane
    end
    
    methods
        function obj = FitPSF_ML_reparam2(psf, par, model)
            if nargin > 2
                % If model is provided, set it
                if strcmpi(model, 'gaussian')
                    obj.model = 'gaussian';
                elseif strcmpi(model, 'hinterer')
                    obj.model = 'hinterer';
                elseif strcmpi(model, 'mortensen')
                    obj.model = 'mortensen';
                end
            end
            
            if nargin > 1
                obj = setInputParameters('FitPSF_ML_reparam2', obj, par);
            end
            
            if nargin > 0
                obj.psf = psf;
                obj.image = psf.image;
                obj.nPhotonEstimate = round(sum(sum(obj.image - obj.noiseEstimate)));
                obj.estimatesPositionDefocus = fitting(obj, model);
            end
        end
        
        %% Fit
        function estimatesPositionDefocus = fitting(obj, model)
            parPsfEstimate = FitPSF_ML_reparam2.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(0, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift;

            % psfEstimate = PSF(parPsfEstimate);
            if strcmpi(obj.model, 'hinterer')
                psfEstimate = PSF(parPsfEstimate);
            elseif strcmpi(obj.model, 'mortensen')
                psfEstimate = PSF_mortensen(parPsfEstimate);
            elseif strcmpi(obj.model, 'gaussian')
                psfEstimate = PSF_gaussian(parPsfEstimate);
            else
                error('Unknown model type: %s', obj.model);
            end

            if strcmpi(model, 'gaussian')
                psfImage = obj.image ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
            elseif strcmpi(model, 'hinterer')
                psfImage = obj.image ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
            elseif strcmpi(model, 'mortensen')
                psfImage = obj.image;% ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
            end

            % estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate, model);%, estimatesPositionDefocus.LS);

        end

        % function estimatesPositionDefocusLS = fitLeastSquaresPSF(obj, image, psfEstimate)
        %     % Define the negative log-likelihood function
        %     funPsf = @(lateralPositionAndDefocus,xdata) createFitPSF(obj, psfEstimate, lateralPositionAndDefocus);
        %     xdata = zeros(obj.psf.nPixels, obj.psf.nPixels);
        %     options = optimoptions('lsqcurvefit','Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 5e-7, 'Display','off');
        % 
        %     startValues = [
        %         obj.parameterStartValues.x.inNanometer, ...
        %         obj.parameterStartValues.y.inNanometer, ...
        %         obj.parameterStartValues.defocus.inNanometer, ...
        %         obj.parameterStartValues.newangle1, ...
        %         obj.parameterStartValues.newangle2, ... 
        %         obj.parameterStartValues.newangle3, ... 
        %         obj.parameterStartValues.photons, ... 
        %         ];
        % 
        %     % Define parameter bounds
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     newangle1Bounds = obj.parameterBounds.newangle1;
        %     newangle2Bounds = obj.parameterBounds.newangle1;
        %     newangle3Bounds = obj.parameterBounds.newangle3;
        %     photonsBounds = obj.parameterBounds.photons;
        % 
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), newangle1Bounds(1), newangle2Bounds(1), newangle3Bounds(1), photonsBounds(1)];
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), newangle1Bounds(2), newangle2Bounds(2), newangle3Bounds(2), photonsBounds(2)];
        % 
        %     estimatesPositionDefocusLS = lsqcurvefit(funPsf, startValues, xdata, image, lowerBounds, upperBounds, options);
        % 
        % end

        function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate, model)
            % Define the negative log-likelihood function
            lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
            costFunction = @(x) -lnpdf(image, x); % Negate likelihood for minimization
        

            if strcmpi(model, 'gaussian') % don't optimise angles in this case

                startValues = [
                    obj.parameterStartValues.x.inNanometer, ...
                    obj.parameterStartValues.y.inNanometer, ...
                    obj.parameterStartValues.defocus.inNanometer, ...
                    obj.parameterStartValues.photons, ... 
                    ];
            
                % Define parameter bounds
                xBounds = obj.parameterBounds.x.inNanometer;
                yBounds = obj.parameterBounds.y.inNanometer;
                defocusBounds = obj.parameterBounds.defocus.inNanometer;
                photonsBounds = obj.parameterBounds.photons;
            
                lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), photonsBounds(1)];
                upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), photonsBounds(2)];
            
                options = optimoptions(@fmincon, 'Display', 'off');
                
                % First set of attempts (3 random attempts)
                fprintf('   Attempt 1\n');
    
                % First optimization attempt with starting values
                firstAttempt = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, [], options);
                firstCost = costFunction(firstAttempt);
                
                % Update best
                bestAttempt = firstAttempt;
                bestCost = firstCost;
                
                % Run 2 more random attempts
                for attempt = 2:3
                    fprintf('   Attempt %d\n', attempt);
                    
                    % Generate random start values
                    newStartValues = startValues;
                    for j = 1:4
                        newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
                    end
                    
                    % Run optimization
                    currentAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, [], options);
                    currentCost = costFunction(currentAttempt);
                    
                    % Update best if this is better
                    if currentCost < bestCost
                        bestAttempt = currentAttempt;
                        bestCost = currentCost;
                    end
                end
                
                % Return the best result and add the cost value
                estimatesPositionDefocusML = bestAttempt;
                estimatesPositionDefocusML(end+1) = bestCost;

            else % optimise angles too for hinterer and mortensen

                % Constraint to ensure cos(az)^2 + sin(az)^2 = 1
                nonlcon = @(x) deal([], x(4)^2 + x(5)^2 + x(6)^2 - 1);
            
                startValues = [
                    obj.parameterStartValues.x.inNanometer, ...
                    obj.parameterStartValues.y.inNanometer, ...
                    obj.parameterStartValues.defocus.inNanometer, ...
                    obj.parameterStartValues.newangle1, ...
                    obj.parameterStartValues.newangle2, ... 
                    obj.parameterStartValues.newangle3, ... 
                    obj.parameterStartValues.photons, ... 
                    ];
            
                % Define parameter bounds
                xBounds = obj.parameterBounds.x.inNanometer;
                yBounds = obj.parameterBounds.y.inNanometer;
                defocusBounds = obj.parameterBounds.defocus.inNanometer;
                newangle1Bounds = obj.parameterBounds.newangle1;
                newangle2Bounds = obj.parameterBounds.newangle2;
                newangle3Bounds = obj.parameterBounds.newangle3;
                photonsBounds = obj.parameterBounds.photons;
            
                lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), newangle1Bounds(1), newangle2Bounds(1), newangle3Bounds(1), photonsBounds(1)];
                upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), newangle1Bounds(2), newangle2Bounds(2), newangle3Bounds(2), photonsBounds(2)];
            
                options = optimoptions(@fmincon, 'Display', 'off');
                
                % First set of attempts (3 random attempts)
                fprintf('   Attempt 1\n');
    
                % First optimization attempt with starting values
                firstAttempt = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
                firstCost = costFunction(firstAttempt);
                firstThetaEst = acos(firstAttempt(6));
                firstThetaEst = mod(firstThetaEst, pi/2);
                
                % Update best
                bestAttempt = firstAttempt;
                bestCost = firstCost;
                bestTheta = firstThetaEst;
                
                % Run 2 more random attempts
                for attempt = 2:3
                    fprintf('   Attempt %d\n', attempt);
                    
                    % Generate random start values
                    newStartValues = startValues;
                    for j = 1:7
                        newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
                    end
                    
                    % Run optimization
                    currentAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
                    currentCost = costFunction(currentAttempt);
                    currentThetaEst = acos(currentAttempt(6));
                    currentThetaEst = mod(currentThetaEst, pi/2);
                    
                    % Update best if this is better
                    if currentCost < bestCost
                        bestAttempt = currentAttempt;
                        bestCost = currentCost;
                        bestTheta = currentThetaEst;
                    end
                end
                
                % Check if bestTheta is within 1 degree of 90 degrees
                if abs(bestTheta - pi/2) < 1*(pi/180)
                    
                    % Run 3 more random attempts
                    for attempt = 4:6
                        fprintf('   Attempt %d\n', attempt);
                        
                        % Generate random start values
                        newStartValues = startValues;
                        for j = 1:7
                            newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
                        end
                        
                        % Run optimization
                        currentAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
                        currentCost = costFunction(currentAttempt);
                        currentThetaEst = acos(currentAttempt(6));
                        currentThetaEst = mod(currentThetaEst, pi/2);
                        
                        % Update best if this is better
                        if currentCost < bestCost
                            bestAttempt = currentAttempt;
                            bestCost = currentCost;
                            bestTheta = currentThetaEst;
                        end
                    end
                    
                    % Check if bestTheta is within 0.001 degrees of 90 degrees
                    if abs(bestTheta - pi/2) < 0.001*(pi/180)
                        
                        % Run up to 5 more attempts or until condition is met
                        maxAttempts = 11; % We've already done 6 attempts (1-6)
                        attempt = 7;
                        
                        while attempt <= maxAttempts && abs(bestTheta - pi/2) < 0.001*(pi/180)
                            fprintf('   Attempt %d\n', attempt);
                            
                            % Generate random start values
                            newStartValues = startValues;
                            for j = 1:7
                                newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
                            end
                            
                            % Run optimization
                            currentAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
                            currentCost = costFunction(currentAttempt);
                            currentThetaEst = acos(currentAttempt(6));
                            currentThetaEst = mod(currentThetaEst, pi/2);
                            
                            % Update best if this is better
                            if currentCost < bestCost
                                bestAttempt = currentAttempt;
                                bestCost = currentCost;
                                bestTheta = currentThetaEst;
                            end
                            
                            attempt = attempt + 1;
                        end
                        
                        % Check if bestTheta > 90 degrees
                        if abs(bestTheta - pi/2) < 0.000001*(pi/180)
                            
                            % Run up to 100 more attempts or until condition is met
                            maxAttempts = attempt + 100; 
                            
                            while attempt <= maxAttempts && bestTheta > pi/2
                                fprintf('   Attempt %d\n', attempt);
                                
                                % Generate random start values
                                newStartValues = startValues;
                                for j = 1:7
                                    newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
                                end
                                
                                % Run optimization
                                currentAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
                                currentCost = costFunction(currentAttempt);
                                currentThetaEst = acos(currentAttempt(6));
                                currentThetaEst = mod(currentThetaEst, pi/2);
                                
                                % Update best if this is better
                                if currentCost < bestCost
                                    bestAttempt = currentAttempt;
                                    bestCost = currentCost;
                                    bestTheta = currentThetaEst;
                                end
                                
                                attempt = attempt + 1;
                            end
                        end
                    end
                end
                
                % Return the best result and add the cost value
                estimatesPositionDefocusML = bestAttempt;
                estimatesPositionDefocusML(end+1) = bestCost;

            end

        end

        function currentlnpdf = lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus) 
            currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
            %currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - log(gamma(z+1)) , 'all');
            % dave feb 2025 - replacing with an all-in-one log gamma function just in case faster
            currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - gammaln(z+1) , 'all');
        end

        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)

            if strcmpi(obj.model, 'gaussian')

                xEstimate = lateralPositionAndDefocus(1);
                yEstimate = lateralPositionAndDefocus(2);
                photonEstimate = lateralPositionAndDefocus(4);
    
                psfEstimate.position = Length([xEstimate, yEstimate, 0], 'nm');
                psfEstimate.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');

            else

                xEstimate = lateralPositionAndDefocus(1);
                yEstimate = lateralPositionAndDefocus(2);
                inclination = acos(lateralPositionAndDefocus(6));
                azimuth = atan2(lateralPositionAndDefocus(5), lateralPositionAndDefocus(4));
                inclination = mod(inclination, pi/2);
                azimuth = mod(azimuth, 2*pi);
                photonEstimate = lateralPositionAndDefocus(7);
    
                psfEstimate.position = Length([xEstimate, yEstimate, 0], 'nm');
                psfEstimate.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');
                psfEstimate.dipole = Dipole(inclination, azimuth); % dave jan 2025 - adding angle optimiser

            end

            % currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
            % for k=1:size(psfEstimate.stageDrift.motion,1)
            %     aberrationCoeffs = getAberrations(psfEstimate,k);
            %     fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
            %     currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
            % end

            % ----------
            % dave jan 2025
            % doing more than the reduced form they were doing            
            if strcmpi(obj.model, 'mortensen')

                % Use Mortensen model

                % Call Python stuff
                % pyDir = pwd; % or specify the exact path where mortensen_simulator.py is located
                pyDir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulation-and-fit/@FitPSF_ML_reparam2'; % or specify the exact path where mortensen_simulator.py is located

                if count(py.sys.path(), pyDir) == 0
                    py.sys.path().insert(int32(0), pyDir);
                end

                % disp(['Sending photon estimate to Python: ', num2str(photonEstimate)]);

                % Run the function defined in your python file
                % currentPsf = py.mortensen_simulator.run_simulator( ...
                currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
                    xEstimate, ... % x
                    yEstimate, ... % y
                    inclination, ... % theta
                    azimuth, ... % phi
                    psfEstimate.nPixels,... % image_size_px
                    double(psfEstimate.pixelSize.inNanometer), ... % pixel_size_nm
                    double(psfEstimate.wavelength.inNanometer), ... % wavelength
                    psfEstimate.refractiveIndices(2), ... % n_objective
                    psfEstimate.refractiveIndices(1), ... % n_sample
                    psfEstimate.objectiveNA, ... % NA
                    photonEstimate ... % n_photons
                );

                % % Show result
                % disp(currentPsf)

                % Convert to matlab array
                % currentPsf = double(currentPsf);
                % Get the shape of the array
                py_shape = currentPsf.shape;
                rows = double(py_shape{1});
                cols = double(py_shape{2});
                
                % Initialize MATLAB array
                psf_matlab = zeros(rows, cols);
                
                % Copy values individually using item() method with explicit integer conversion
                for i = 0:(rows-1)
                    for j = 0:(cols-1)
                        % Use py.int to explicitly convert indices to Python integers
                        psf_matlab(i+1, j+1) = double(currentPsf.item(py.tuple({py.int(i), py.int(j)})));
                    end
                end
                
                currentPsf = psf_matlab;




                % this bit is the equivalent of the getIntensitiesCamera()
                % stuff done in hinterer
                totalIntensity = sum(sum(currentPsf));
                currentPsf = currentPsf / totalIntensity * photonEstimate;
                currentFitPSF = currentPsf;

            else

                % Select appropriate BackFocalPlane function based on model
                if strcmpi(obj.model, 'gaussian')
                    bfp = BackFocalPlane_gaussian(psfEstimate); % Use Gaussian model
                elseif strcmpi(obj.model, 'hinterer')
                    bfp = BackFocalPlane(psfEstimate); % Use Hinterer model
                end
    
                    psfEstimate.backFocalPlane = bfp;
        
                    % dave
                    psfEstimate.fieldBFP.x = bfp.electricField.x;
                    psfEstimate.fieldBFP.y = bfp.electricField.y;
        
                    % % Apply phase mask
                    % psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
                    % psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);
        
                    % % Apply attenuation mask
                    % psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
                    % psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);
        
                    currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
                    for k=1:size(psfEstimate.stageDrift.motion,1)
                        % Apply aberrations
                        aberrationCoeffs = getAberrations(psfEstimate,k);
                        fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
                        % Get image from BFP field
                        % currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
                        currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
                    end

                    totalIntensity = sum(currentPsf,'all');
                    currentPsf = currentPsf ./ totalIntensity * photonEstimate + obj.noiseEstimate;
                    % disp(currentPsf)
                    currentFitPSF = currentPsf ./ norm(currentPsf);
                    % disp(currentFitPSF)

            end
        
            currentFitPSF = adjustExcitation(psfEstimate, currentFitPSF);
            currentFitPSF = applyShotNoise(psfEstimate, currentFitPSF);
            currentFitPSF = addBackgroundNoise(psfEstimate, currentFitPSF);

            % % Output as tif
            % counter = 1000*rand();
            % output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_45/fitting/sim_frame%06d.tif', round(counter));
            % psf_total_image = uint32(currentFitPSF);
            % t = Tiff(output_path, 'w');
            % tagstruct.ImageLength = 19;  % Set image height
            % tagstruct.ImageWidth = 19;   % Set image width
            % tagstruct.Photometric = Tiff.Photometric.MinIsBlack;  % Grayscale
            % tagstruct.BitsPerSample = 32;  % 16-bit per pixel (or 32-bit if needed)
            % tagstruct.SamplesPerPixel = 1;  % Grayscale (1 channel)
            % tagstruct.RowsPerStrip = 16;   % Strip length for compression
            % tagstruct.Compression = Tiff.Compression.LZW;  % Lossless compression (optional)
            % tagstruct.Software = 'MATLAB';
            % t.setTag(tagstruct);
            % t.write(psf_total_image);
            % t.close();

        end

    end

    methods (Static)
        par = readParametersEstimate(psf);
    end
end