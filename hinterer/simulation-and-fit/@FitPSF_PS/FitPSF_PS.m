classdef FitPSF_PS
    
    properties
        psf PSF 
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
            'inclination', [0, pi/2], ...       % dave jan 2025 - adding angle optimiser
            'azimuth', [-0.0174533, (2*pi)+0.0174533]);          % dave jan 2025 - adding angle optimiser
        
        parameterStartValues = struct( ...
            'x', Length(-100 + 200 * rand(), 'nm'), ...
            'y', Length(-100 + 200 * rand(), 'nm'), ...
            'defocus', Length(-500 + 1000 * rand(), 'nm'), ...
            'inclination', rand() * pi/2, ...   % dave jan 2025 - adding angle optimiser
            'azimuth', rand() * 2 * pi ...  % dave jan 2025 - adding angle optimiser
            );
        
        % Fit result
        estimatesPositionDefocus
    end
    
    methods
        function obj = FitPSF_PS(psf, par)
            if nargin > 1
                obj = setInputParameters('FitPSF_PS', obj, par);
            end
            if nargin > 0
                obj.psf = psf;
                obj.image = psf.image;
                obj.nPhotonEstimate = round(sum(sum(obj.image - obj.noiseEstimate)));
                obj.estimatesPositionDefocus = fitting(obj);
            end
        end
        
        %% Fit
        function estimatesPositionDefocus = fitting(obj)
            parPsfEstimate = FitPSF_PS.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(0, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift; 

            psfEstimate = PSF(parPsfEstimate);

            psfImage = obj.image ./ norm(obj.image);
            
            % plotObjectiveFunction(obj, psfEstimate, psfImage);

            % estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);
            % estimatesPositionDefocus.GA = fitGeneticAlgorithmPSF(obj, psfImage, psfEstimate); % dave jan 2025
            % estimatesPositionDefocus.SA = fitSimulatedAnnealingPSF(obj, psfImage, psfEstimate); % dave jan 2025
            estimatesPositionDefocus.PS = fitParticleSwarmPSF(obj, psfImage, psfEstimate); % dave jan 2025

            % % dave jan 2025
            % % hybrid optimisation: ML for position, GA for angles
            % estimatesPositionDefocus.LS = fitLeastSquaresPSF_noinc(obj, psfImage, psfEstimate);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF_noinc(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);
            % estimatesPositionDefocus.GA = fitGeneticAlgorithmPSF_onlyinc(obj, psfImage, psfEstimate, estimatesPositionDefocus.ML);
            % estimatesPositionDefocus.SA = fitSimulatedAnnealingPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.ML);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF_noinc(obj, psfImage, psfEstimate, estimatesPositionDefocus.GA);
        end
        
        % function estimatesPositionDefocusLS = fitLeastSquaresPSF(obj, image, psfEstimate)
        %     funPsf = @(lateralPositionAndDefocus,xdata) createFitPSF(obj, psfEstimate, lateralPositionAndDefocus);
        %     xdata = zeros(obj.psf.nPixels,obj.psf.nPixels);
        %     options = optimoptions('lsqcurvefit','Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 5e-7, 'Display','off');
        % 
        %     startValues = [obj.parameterStartValues.x.inNanometer, ...
        %         obj.parameterStartValues.y.inNanometer, ...
        %         obj.parameterStartValues.defocus.inNanometer, ...
        %         obj.parameterStartValues.azimuth]; % dave jan 2025 - adding angle optimiser
        % 
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     azimuthBounds = obj.parameterBounds.azimuth;     % dave jan 2025 - adding angle optimiser
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)]; % dave jan 2025 - adding angle optimiser, unfixed az
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)]; % dave jan 2025 - adding angle optimiser, unfixed az
        %     estimatesPositionDefocusLS = lsqcurvefit(funPsf, startValues, xdata, image, lowerBounds, upperBounds, options);
        % end

        % function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate, startValues)
        %     lnpdf = @(z,lateralPositionAndDefocus) lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus);
        %     options = optimoptions(@fmincon, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
        %     % options = optimoptions(@fmincon, ...
        %     %     'Display', 'off', ...               % Do not display output
        %     %     'StepTolerance', 1e-20, ...          % Stop when the step size is less than 1e-6
        %     %     'OptimalityTolerance', 1e-20, ...    % Stop when the gradient is less than 1e-6
        %     %     'MaxIterations', 10000000, ...          % Stop after 1000 iterations
        %     %     'MaxFunctionEvaluations', 50000000, ... % Stop after 5000 function evaluations
        %     %     'FunctionTolerance', 1e-20);        % Stop if the function value change is less than 1e-6
        %     % estimatesPositionDefocusML = fminunc(@(x) -lnpdf(image, x), startValues, options);
        %     % dave jan 2025
        %     % using constrained version, because why wouldnt you?
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     azimuthBounds = obj.parameterBounds.azimuth;
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)]; % unfixed azimuth
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)]; % unfixed azimuth
        %     estimatesPositionDefocusML = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, [], options);
        % end

        % % dave jan 2025 - trying better optimiser to avoid local minima
        % function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate, startValues)
        %     lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
        %     options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10); % Using optimset for fminsearch
        %     estimatesPositionDefocusML = fminsearch(@(x) -lnpdf(image, x), startValues, options);
        % end

        function estimatesPositionDefocusPS = fitParticleSwarmPSF(obj, image, psfEstimate)
            % Define the negative log-likelihood function
            lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
            costFunction = @(x) -lnpdf(image, x); % Negate likelihood for minimization
        
            % Define parameter bounds
            defocusBounds = obj.parameterBounds.defocus.inNanometer;
            xBounds = obj.parameterBounds.x.inNanometer;
            yBounds = obj.parameterBounds.y.inNanometer;
            inclinationBounds = obj.parameterBounds.inclination;
            azimuthBounds = obj.parameterBounds.azimuth;

            lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), inclinationBounds(1), azimuthBounds(1)];
            upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), inclinationBounds(2), azimuthBounds(2)];

            % Set options for the Particle Swarm Optimization

            % Create a uniform grid for initial swarm positions
            % numParticlesPerDim = 3;
            [xGrid, yGrid, defocusGrid, inclinationGrid, azimuthGrid] = ndgrid(...
                linspace(xBounds(1), xBounds(2), 3), ...
                linspace(yBounds(1), yBounds(2), 3), ...
                linspace(defocusBounds(1), defocusBounds(2), 2), ...
                linspace(inclinationBounds(1), inclinationBounds(2), 2), ...
                linspace(azimuthBounds(1), azimuthBounds(2), 2));
            % Reshape grid into an initial swarm matrix
            initialSwarm = [xGrid(:), yGrid(:), defocusGrid(:), inclinationGrid(:),azimuthGrid(:)];
            numParticles = 72; % this needs to be the multiplicative sum of each param given in the linspace above
            initialSwarm = initialSwarm(1:numParticles, :);
            
            options = optimoptions('particleswarm', ...
                'Display', 'off', ...         % Show optimization progress
                'SwarmSize', numParticles, ...          % Number of particles
                'MaxIterations', 1e20, ...      % Maximum number of iterations
                'InitialSwarmMatrix', initialSwarm, ...      % Swarm distribution
                'FunctionTolerance', 1e-9);   % Stopping criteria for function improvement

            % Run Particle Swarm Optimization
            numVariables = numel(lowerBounds);
            estimatesPositionDefocusPS = particleswarm(costFunction, numVariables, lowerBounds, upperBounds, options);

            % add in something to save the value of the objetive function
            estimatesPositionDefocusPS(end+1) = costFunction(estimatesPositionDefocusPS);
        end



        function currentlnpdf = lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus) 
            currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
            currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - log(gamma(z+1)) , 'all');
        end

        
        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)
            psfEstimate.position = Length([lateralPositionAndDefocus(1:2), 0], 'nm');
            psfEstimate.defocus = Length(lateralPositionAndDefocus(3), 'nm');
            % lateralPositionAndDefocus(4) = mod(lateralPositionAndDefocus(4), pi/2);
            % lateralPositionAndDefocus(5) = mod(lateralPositionAndDefocus(5), 2*pi);

            psfEstimate.dipole = Dipole(lateralPositionAndDefocus(4), lateralPositionAndDefocus(5)); % dave jan 2025 - adding angle optimiser

            % currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
            % for k=1:size(psfEstimate.stageDrift.motion,1)
            %     aberrationCoeffs = getAberrations(psfEstimate,k);
            %     fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
            %     currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
            % end
            % totalIntensity = sum(currentPsf,'all');
            % currentPsf = currentPsf ./ totalIntensity * obj.nPhotonEstimate + obj.noiseEstimate;
            % currentFitPSF = currentPsf ./ norm(currentPsf);

            % dave jan 2025
            % doing more than the reduced form they were doing            
            
            bfp = BackFocalPlane(psfEstimate);
            % bfp = BackFocalPlane_gaussian(psfEstimate); % use this if want just Gaussian
            psfEstimate.backFocalPlane = bfp;

            % Apply phase mask
            psfEstimate.fieldBFP.x = psfEstimate.phaseMaskObj.apply(bfp.electricField.x);
            psfEstimate.fieldBFP.y = psfEstimate.phaseMaskObj.apply(bfp.electricField.y);

            % Apply attenuation mask
            psfEstimate.fieldBFP.x = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.x);
            psfEstimate.fieldBFP.y = psfEstimate.attenuationMaskObj.apply(psfEstimate.fieldBFP.y);

            currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
            for k=1:size(psfEstimate.stageDrift.motion,1)
                % Apply aberrations
                aberrations = getAberrations(psfEstimate,k);
                aberratedFieldBFP = applyAberrations(psfEstimate, aberrations);
                
                % Get image from BFP field
                currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, aberratedFieldBFP)./size(psfEstimate.stageDrift.motion,1);
            end

            currentPsf = adjustExcitation(psfEstimate, currentPsf);
            currentPsf = applyShotNoise(psfEstimate, currentPsf);
            currentPsf = addBackgroundNoise(psfEstimate, currentPsf);

            totalIntensity = sum(currentPsf,'all');
            currentPsf = currentPsf ./ totalIntensity * obj.nPhotonEstimate + obj.noiseEstimate;
            currentFitPSF = currentPsf ./ norm(currentPsf);

            % dave jan 2025
            % disp(['X: ', num2str(lateralPositionAndDefocus(1))]);
            % disp(['Inclination: ', num2str(lateralPositionAndDefocus(4))]);

            % % dave jan 2025 - print to check if PSF is updated every time
            % persistent iterationCounter;  % Keeps track of iterations
            % if isempty(iterationCounter)
            %     iterationCounter = 0;
            % end
            % iterationCounter = iterationCounter + 1;  % Increment counter
            % if mod(iterationCounter, 10) == 0  % Only output on every 10th iteration
            %     outputDirectory = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/optimiser_output';  % Define the output directory
            %     if ~exist(outputDirectory, 'dir')
            %         mkdir(outputDirectory);  % Create the directory if it doesn't exist
            %     end
            %     timestamp = datestr(now, 'yyyymmdd_HHMMSS_FFF'); 
            %     filename = sprintf('iteration_%s.png', timestamp);
            %     imwrite(mat2gray(currentFitPSF), fullfile(outputDirectory, filename));
            % end

        end

        % heat map objective function plot
        function plotObjectiveFunction(obj, psfEstimate, image)
            % Define radial and azimuthal range
            rVals = linspace(0, 60, 60);           % Radial values (distance from origin)
            azimuthVals = linspace(0, 2*pi, 180);  % Azimuthal values (angle)
        
            % Fix defocus to 0
            defocus = 0;
        
            % Initialize cost function values
            costMatrix = zeros(length(rVals), length(azimuthVals));
        
            % Compute cost function over radial and azimuthal coordinates
            for i = 1:length(rVals)
                for j = 1:length(azimuthVals)
                    % Convert (r, θ) to (x, y)
                    x = rVals(i) * cos(azimuthVals(j));
                    y = rVals(i) * sin(azimuthVals(j));
                    params = [x, y, defocus, azimuthVals(j)];
                    costMatrix(i, j) = -obj.lnpdfFunction(psfEstimate, image, params); % Negative log-likelihood
                end
            end
        
            % Convert polar to Cartesian for 3D surface plotting
            [AZI, R] = meshgrid(azimuthVals, rVals);
            [X, Y] = pol2cart(AZI, R);  % Convert polar to Cartesian coordinates
        
            % Create a 3D surface plot
            figure;
            surf(X, Y, costMatrix, 'EdgeColor', 'none'); % Remove edges for smooth look
            view(3); % 3D view
            colorbar;
            colormap jet;
            xlabel('X Position (nm)');
            ylabel('Y Position (nm)');
            zlabel('Objective Function Value');
            title('Objective Function in Polar Coordinates (Defocus = 0)');
        end
        
        %% Plot
        % function plot(obj)
        %     plot(obj.psf)
        %     hold on
        %     size = 12;
        %     width = 2.5;
        %     center = (obj.psf.nPixels+1)/2;
        %     plot(center+obj.estimatesPositionDefocus.LS(1)/100, center+obj.estimatesPositionDefocus.LS(2)/100,'Marker','o','MarkerSize',size,'Color','black','LineWidth', width)
        %     plot(center+obj.estimatesPositionDefocus.ML(1)/100, center+obj.estimatesPositionDefocus.ML(2)/100,'Marker','+','MarkerSize',size,'Color',[1 1 1]*0.8,'LineWidth', width)
        %     plot(obj.psf.positionInPixelFromOrigin(1), obj.psf.positionInPixelFromOrigin(2),'Marker','x','MarkerSize',size,'Color','red','LineWidth', width)
        %     axis equal
        %     axis tight
        %     cb = colorbar;
        %     ylabel(cb,'Intensity','FontSize',15)
        % end
    end

    methods (Static)
        par = readParametersEstimate(psf);
    end
end

