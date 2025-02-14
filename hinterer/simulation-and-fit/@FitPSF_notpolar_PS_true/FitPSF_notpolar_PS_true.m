classdef FitPSF_notpolar_PS_true
    
    properties
        psf PSF 
        image
        angleInclinationEstimate (1,1) {mustBeInFullRadialRange} % dave jan 2025 - commented out for adding angle optimiser
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
        function obj = FitPSF_notpolar_PS_true(psf, par)
            if nargin > 1
                obj = setInputParameters('FitPSF_notpolar_PS_true', obj, par);
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
            parPsfEstimate = FitPSF_notpolar_PS_true.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(obj.parameterStartValues.inclination, obj.parameterStartValues.azimuth);
            parPsfEstimate.position = Length([obj.parameterStartValues.x.inNanometer obj.parameterStartValues.y.inNanometer obj.parameterStartValues.defocus.inNanometer], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0.0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift; 

            psfEstimate = PSF(parPsfEstimate);

            psfImage = obj.image ./ norm(obj.image);
            
            % plotObjectiveFunction(obj, psfEstimate, psfImage);

            % estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);
            % estimatesPositionDefocus.GA = fitGeneticAlgorithmPSF(obj, psfImage, psfEstimate); % dave jan 2025
            % estimatesPositionDefocus.SA = fitSimulatedAnnealingPSF(obj, psfImage, psfEstimate); % dave jan 2025
            estimatesPositionDefocus.obj_func_true = fitParticleSwarmPSF(obj, psfImage, psfEstimate); % dave jan 2025

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

        % just evaluate objective function at ground truth
        function estimatesPositionDefocusPS = fitParticleSwarmPSF(obj, image, psfEstimate)

            groundTruthParams = [obj.parameterStartValues.x.inNanometer, ...
                obj.parameterStartValues.y.inNanometer, ...
                obj.parameterStartValues.defocus.inNanometer, ...
                obj.parameterStartValues.inclination, ...
                obj.parameterStartValues.azimuth];

            currentPSF = createFitPSF(obj, psfEstimate, groundTruthParams);
            costValue = -sum(image .* log(currentPSF) - currentPSF - log(gamma(image + 1)), 'all');

            estimatesPositionDefocusPS = costValue;
        end
        
        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)
            psfEstimate.position = Length([lateralPositionAndDefocus(1:2), 0], 'nm');
            psfEstimate.defocus = Length(lateralPositionAndDefocus(3), 'nm');
            lateralPositionAndDefocus(4) = mod(lateralPositionAndDefocus(4), 2*pi);

            psfEstimate.dipole = Dipole(obj.angleInclinationEstimate, lateralPositionAndDefocus(4)); % dave jan 2025 - adding angle optimiser

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

        % % heat map objective function plot
        % function plotObjectiveFunction(obj, psfEstimate, image)
        % 
        %     tic;
        % 
        %     % Define radial and azimuthal range
        %     rVals = linspace(0, 60, 120);           % Radial values (distance from origin)
        %     azimuthVals = linspace(0, 2*pi, 360);  % Azimuthal values (angle)
        % 
        %     % Fix defocus to 0
        %     defocus = 0;
        % 
        %     % Initialize cost function values
        %     costMatrix = zeros(length(rVals), length(azimuthVals));
        % 
        %     % Compute cost function over radial and azimuthal coordinates
        %     for i = 1:length(rVals)
        %         for j = 1:length(azimuthVals)
        %             % tic;
        %             % Convert (r, θ) to (x, y)
        %             x = rVals(i) * cos(azimuthVals(j));
        %             y = rVals(i) * sin(azimuthVals(j));
        %             params = [x, y, defocus, azimuthVals(j)];
        %             currentPSF = createFitPSF(obj, psfEstimate, params);
        %             costValue = -sum(image .* log(currentPSF) - currentPSF - log(gamma(image + 1)), 'all');
        %             costMatrix(i, j) = costValue; % Negative log-likelihood
        %             % elapsed_time = toc;
        %             % fprintf('    Time to plot: %.20f seconds\n', elapsed_time);
        %         end
        %     end
        % 
        %     % Convert polar to Cartesian for 3D surface plotting
        %     [AZI, R] = meshgrid(azimuthVals, rVals);
        %     [X, Y] = pol2cart(AZI, R);  % Convert polar to Cartesian coordinates
        % 
        %     % Create a 3D surface plot
        %     figure;
        %     surf(X, Y, costMatrix, 'EdgeColor', 'none'); % Remove edges for smooth look
        %     view(3); % 3D view
        %     colorbar;
        %     colormap jet;
        %     xlabel('x, nm)');
        %     ylabel('y, nm)');
        %     zlabel('obj func');
        %     title('Objective Func (θ=0)');
        % 
        %     elapsed_time = toc;
        %     fprintf('    Time to plot: %.2f seconds\n', elapsed_time);
        % 
        % end

        function plotObjectiveFunction(obj, psfEstimate, image)
                
            xVals = linspace(-50, 50, 50);
            yVals = linspace(-50, 50, 50);
            
            % Fix defocus and azimuth to 0
            defocus = obj.parameterStartValues.defocus.inNanometer;
            inclination = obj.parameterStartValues.inclination;
            azimuth = obj.parameterStartValues.azimuth;

            % Initialize cost function values
            costMatrix = zeros(length(xVals), length(yVals));

            % Compute cost function for each inclination value and place on subplot
        
            tic;

            % Initialize cost function values
            costMatrix = zeros(length(xVals), length(yVals));
    
            % Compute the cost function over x and y coordinates
            for i = 1:length(xVals)
                for j = 1:length(yVals)
                    % Current (x, y) position
                    x = xVals(i);
                    y = yVals(j);
                    params = [x, y, defocus, inclination, azimuth];  % Azimuth is set to 0
                    currentPSF = createFitPSF(obj, psfEstimate, params);
                    costValue = -sum(image .* log(currentPSF) - currentPSF - log(gamma(image + 1)), 'all');
                    costMatrix(i, j) = costValue;  % Negative log-likelihood
                end
            end
    
            % Create a 3D surface plot
            figure;
            surf(xVals, yVals, costMatrix, 'EdgeColor', 'none');
            view(3);
            colorbar;
            colormap jet;
            xlabel('x, nm)');
            ylabel('y, nm)');
            zlabel('obj func');
            title(sprintf('obj func (θ = %.2f°)', inclination*180/pi));
    
            elapsed_time = toc;
            fprintf('    Time to plot: %.2f seconds\n', elapsed_time);

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

