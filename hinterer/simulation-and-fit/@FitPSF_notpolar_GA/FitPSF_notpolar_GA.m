classdef FitPSF_notpolar_GA
    
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
            'azimuth', [-Inf, Inf]);          % dave jan 2025 - adding angle optimiser
        
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
        function obj = FitPSF_notpolar_GA(psf, par)
            if nargin > 1
                obj = setInputParameters('FitPSF_notpolar_GA', obj, par);
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
            parPsfEstimate = FitPSF_notpolar_GA.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(obj.angleInclinationEstimate, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift; 

            psfEstimate = PSF(parPsfEstimate);

            psfImage = obj.image ./ norm(obj.image);
            
            % estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);
            estimatesPositionDefocus.GA = fitGeneticAlgorithmPSF(obj, psfImage, psfEstimate); % dave jan 2025
            % estimatesPositionDefocus.SA = fitSimulatedAnnealingPSF(obj, psfImage, psfEstimate); % dave jan 2025

            % % dave jan 2025
            % % hybrid optimisation: ML for position, GA for angles
            % estimatesPositionDefocus.LS = fitLeastSquaresPSF_noinc(obj, psfImage, psfEstimate);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF_noinc(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);
            % estimatesPositionDefocus.GA = fitGeneticAlgorithmPSF_onlyinc(obj, psfImage, psfEstimate, estimatesPositionDefocus.ML);
            % estimatesPositionDefocus.SA = fitSimulatedAnnealingPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.ML);
            % estimatesPositionDefocus.ML = fitMaxLikelihoodPSF_noinc(obj, psfImage, psfEstimate, estimatesPositionDefocus.GA);
        end
        
        function estimatesPositionDefocusLS = fitLeastSquaresPSF(obj, image, psfEstimate)
            funPsf = @(lateralPositionAndDefocus,xdata) createFitPSF(obj, psfEstimate, lateralPositionAndDefocus);
            xdata = zeros(obj.psf.nPixels,obj.psf.nPixels);
            options = optimoptions('lsqcurvefit','Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 5e-7, 'Display','off');

            startValues = [obj.parameterStartValues.x.inNanometer, ...
                obj.parameterStartValues.y.inNanometer, ...
                obj.parameterStartValues.defocus.inNanometer, ...
                obj.parameterStartValues.azimuth]; % dave jan 2025 - adding angle optimiser

            defocusBounds = obj.parameterBounds.defocus.inNanometer;
            xBounds = obj.parameterBounds.x.inNanometer;
            yBounds = obj.parameterBounds.y.inNanometer;
            azimuthBounds = obj.parameterBounds.azimuth;     % dave jan 2025 - adding angle optimiser
            lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)]; % dave jan 2025 - adding angle optimiser, unfixed az
            upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)]; % dave jan 2025 - adding angle optimiser, unfixed az
            estimatesPositionDefocusLS = lsqcurvefit(funPsf, startValues, xdata, image, lowerBounds, upperBounds, options);
        end

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

        % dave jan 2025 - trying better optimiser to avoid local minima
        % function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate, startValues)
        %     lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
        %     options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10); % Using optimset for fminsearch
        %     estimatesPositionDefocusML = fminsearch(@(x) -lnpdf(image, x), startValues, options);
        % end

        % % dave jan 2025 - no inclination, for hybrid optimising
        % function estimatesPositionDefocusLS = fitLeastSquaresPSF_noinc(obj, image, psfEstimate)
        %     % Fix inclination to its initial value
        %     fixedInclination = obj.parameterStartValues.inclination; % Fix inclination
        % 
        %     % Define the fitting function, excluding inclination
        %     funPsf = @(lateralPositionDefocusAzimuth, xdata) createFitPSF(...
        %         obj, psfEstimate, [lateralPositionDefocusAzimuth(1:3), fixedInclination, lateralPositionDefocusAzimuth(4)]);
        % 
        %     % Dummy xdata for lsqcurvefit
        %     xdata = zeros(obj.psf.nPixels, obj.psf.nPixels);
        % 
        %     % Optimization options
        %     options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', ...
        %         'OptimalityTolerance', 5e-7, 'Display', 'off');
        % 
        %     % Start values excluding inclination
        %     startValues = [obj.parameterStartValues.x.inNanometer, ...
        %                    obj.parameterStartValues.y.inNanometer, ...
        %                    obj.parameterStartValues.defocus.inNanometer, ...
        %                    obj.parameterStartValues.azimuth]; % Exclude inclination
        % 
        %     % Bounds excluding inclination
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     azimuthBounds = obj.parameterBounds.azimuth;
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)];
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)];
        % 
        %     % Perform optimization
        %     reducedEstimates = lsqcurvefit(funPsf, startValues, xdata, image, lowerBounds, upperBounds, options);
        % 
        %     % Combine the fixed inclination with the optimized values
        %     estimatesPositionDefocusLS = [reducedEstimates(1:3), fixedInclination, reducedEstimates(4)];
        % end
        % 
        % % dave jan 2025 - no inclination, for hybrid optimising
        % function estimatesPositionDefocusML = fitMaxLikelihoodPSF_noinc(obj, image, psfEstimate, startValues)
        %     fixedInclination = startValues(4);  % The initial inclination value
        %     lnpdf = @(z,lateralPositionAndDefocus) lnpdfFunction(obj,psfEstimate,z,[lateralPositionAndDefocus(1:3), fixedInclination, lateralPositionAndDefocus(4)]);
        % 
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     azimuthBounds = obj.parameterBounds.azimuth;
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)];
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)];
        % 
        %     reducedStartValues = [startValues(1:3), startValues(5)];
        %     options = optimoptions(@fmincon, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
        %     reducedEstimates = fmincon(@(x) -lnpdf(image, x), reducedStartValues, [], [], [], [], lowerBounds, upperBounds, [], options);
        %     estimatesPositionDefocusML = [reducedEstimates(1:3), fixedInclination, reducedEstimates(4)];
        % end

        % % dave jan 2025 - add inclination penalty to fminunc
        % % no effect
        % function estimatesPositionDefocusMLpenalty = fitMaxLikelihoodPSFpenalty(obj, image, psfEstimate, startValues)
        %     lnpdf = @(z,lateralPositionAndDefocus) lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus);
        %     options = optimoptions(@fminunc, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
        %     % -lnpdf(image,x) changed to:
        %     % penalty_weight*inclination^2 - lnpdf(image,x)
        %     estimatesPositionDefocusML = fminunc(@(x) (1e-40)*log(1 + 1/abs(x(4))) - lnpdf(image,x), startValues, options);
        % end

        % dave jan 2025 - try genetic algorithm instead of fminunc
        % performance of this was pretty bad for all params, 
        % and with longer computation time than fminunc
        % (which was only bad for inclination + azimuth)
        function estimatesPositionDefocusGA = fitGeneticAlgorithmPSF(obj, image, psfEstimate)
            % Define the cost function to minimize
            costFunction = @(x) lnpdfFunction(obj, psfEstimate, image, x);

            % Lower and Upper bounds for each of the parameters        
            defocusBounds = obj.parameterBounds.defocus.inNanometer;
            xBounds = obj.parameterBounds.x.inNanometer;
            yBounds = obj.parameterBounds.y.inNanometer;
            azimuthBounds = obj.parameterBounds.azimuth;     % dave jan 2025 - adding angle optimiser
            lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)]; % dave jan 2025 - adding angle optimiser
            upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)]; % dave jan 2025 - adding angle optimiser


            % Set GA options (make sure they are correctly configured)
            options = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 100, 'MaxGenerations', 15, ...
                                    'EliteCount', 10, 'CrossoverFraction', 0.9, 'PlotFcn', @gaplotbestf);

            % Run the genetic algorithm to minimize the cost function
            [optimizedParams, ~] = ga(costFunction, 4, [], [], [], [], lowerBounds, upperBounds, [], [], options);

            % Store the result in the estimatesPositionDefocusGA structure
            estimatesPositionDefocusGA = optimizedParams;  % Store the parameters directly
        end

        % % dave jan 2025 - try simulated annealing instead of fminunc
        % % performance of this was pretty bad for all params, 
        % % and with longer computation time than fminunc
        % % (which was only bad for inclination + azimuth)
        % function estimatesPositionDefocusSA = fitSimulatedAnnealingPSF(obj, image, psfEstimate)
        %     % Define the cost function to minimize
        %     costFunction = @(x) lnpdfFunction(obj, psfEstimate, image, x);
        % 
        %     startValues = [obj.parameterStartValues.x.inNanometer, ...
        %         obj.parameterStartValues.y.inNanometer, ...
        %         obj.parameterStartValues.defocus.inNanometer, ...
        %         obj.parameterStartValues.inclination, ... % dave jan 2025 - adding angle optimiser
        %         obj.parameterStartValues.azimuth]; % dave jan 2025 - adding angle optimiser
        % 
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     inclinationBounds = obj.parameterBounds.inclination; % dave jan 2025 - adding angle optimiser
        %     azimuthBounds = obj.parameterBounds.azimuth;     % dave jan 2025 - adding angle optimiser
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), inclinationBounds(1), azimuthBounds(1)]; % dave jan 2025 - adding angle optimiser
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), inclinationBounds(2), azimuthBounds(2)]; % dave jan 2025 - adding angle optimiser
        % 
        % 
        %     % Set SA options (make sure they are correctly configured)
        %     % options = optimoptions('simulannealbnd', 'Display', 'iter', 'MaxIterations', 1000, ...
        %     %                        'InitialTemperature', 100, 'TemperatureFcn', @temperatureexp, 'ReannealInterval', 50);
        %     options = optimoptions('simulannealbnd', ...
        %         'Display', 'iter', ...               % Show iteration details in the command window
        %         'MaxIterations', 10000, ...            % Maximum number of iterations (higher for more thorough search)
        %         'MaxFunctionEvaluations', 10000, ...  % Maximum number of function evaluations (use if large function evaluations are expensive)
        %         'InitialTemperature', 10000, ...        % Initial "temperature" to explore the solution space
        %         'TemperatureFcn', @temperatureexp, ...        % Initial "temperature" to explore the solution spaceQ
        %         'FunctionTolerance', 1e-20);
        % 
        %     % Run the simulated annealing to minimize the cost function
        %     [optimizedParams, ~] = simulannealbnd(costFunction, startValues, lowerBounds, upperBounds, options);
        % 
        %     % Store the result in the estimatesPositionDefocusGA structure
        %     estimatesPositionDefocusSA = optimizedParams;  % Store the parameters directly
        % end
        % 
        % dave jan 2025
        % new GA that only goes on angles. ML used for positions
        % function estimatesPositionDefocusGA = fitGeneticAlgorithmPSF_onlyinc(obj, image, psfEstimate, optimizedParams)
        %     % Define the cost function to minimize
        %     costFunction = @(x) lnpdfFunction(obj, psfEstimate, image, x);
        % 
        %     % Lower and Upper bounds for each of the parameters (only inclination and azimuth)
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     inclinationBounds = obj.parameterBounds.inclination; 
        %     azimuthBounds = obj.parameterBounds.azimuth;
        % 
        %     % Fixed parameters from the first optimization (x, y, defocus)
        %     fixedX = optimizedParams(1);
        %     fixedY = optimizedParams(2);
        %     fixedDefocus = optimizedParams(3);
        % 
        %     % The initial guess for the GA is just the initial values of inclination and azimuth
        %     initialGuessInc = optimizedParams(4);  % inclination and azimuth
        %     initialGuessAz = optimizedParams(5);  % inclination and azimuth
        % 
        %     % Lower and upper bounds for inclination and azimuth
        %     lowerBounds = [inclinationBounds(1), azimuthBounds(1)];
        %     upperBounds = [inclinationBounds(2), azimuthBounds(2)];
        % 
        %     % Set GA options (make sure they are correctly configured)
        %     options = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 100, 'MaxGenerations', 20, ...
        %                             'EliteCount', 50, 'CrossoverFraction', 0.9, 'PlotFcn', @gaplotbestf);
        % 
        %     % Run the genetic algorithm to minimize the cost function (the 1 here is number of params to optimise over)
        %     [optimizedInclinationAzimuth, ~] = ga(@(params) costFunction([fixedX, fixedY, fixedDefocus, params(1), params(2)]), 2, [], [], [], [], lowerBounds, upperBounds, [], options);
        % 
        %     % Store the result in the estimatesPositionDefocusGA structure
        %     estimatesPositionDefocusGA = [fixedX, fixedY, fixedDefocus, optimizedInclinationAzimuth(1), optimizedInclinationAzimuth(2)];
        % end
        % 
        %
        % % dave jan 2025
        % % new SA that only goes on angles. ML used for positions
        % function estimatesPositionDefocusSA = fitSimulatedAnnealingPSF(obj, image, psfEstimate, optimizedParams)
        %     % Define the cost function to minimize (same as for the GA)
        %     costFunction = @(x) lnpdfFunction(obj, psfEstimate, image, x);
        % 
        %     % Fixed parameters from the first optimization (x, y, defocus)
        %     fixedX = optimizedParams(1);
        %     fixedY = optimizedParams(2);
        %     fixedDefocus = optimizedParams(3);
        %     fixedAzimuth = optimizedParams(5); % Keep azimuth fixed (for now)
        % 
        %     % The initial guess for the SA is just the initial value of inclination
        %     initialGuess = optimizedParams(4);  % inclination
        % 
        %     % Lower and upper bounds for inclination
        %     inclinationBounds = obj.parameterBounds.inclination; 
        %     lowerBounds = [inclinationBounds(1)];
        %     upperBounds = [inclinationBounds(2)];
        % 
        %     % Define the objective function for inclination optimization (fixed x, y, defocus, azimuth)
        %     objectiveFunction = @(inclination) costFunction([fixedX, fixedY, fixedDefocus, inclination, fixedAzimuth]);
        % 
        %     % Set Simulated Annealing options
        %     % optionsSA = optimoptions('simulannealbnd', 'Display', 'iter', 'MaxIterations', 1000, 'MaxTime', 6000);
        %     optionsSA = optimoptions('simulannealbnd', 'Display', 'iter', 'MaxIterations', 1000, ...
        %         'InitialTemperature', 100, 'TemperatureFcn', @temperatureexp, 'ReannealInterval', 50);
        %     % Run Simulated Annealing to minimize the cost function for inclination
        %     [optimizedInclination, ~] = simulannealbnd(objectiveFunction, initialGuess, lowerBounds, upperBounds, optionsSA);
        % 
        %     % Combine the optimized inclination with the fixed x, y, defocus, and azimuth
        %     estimatesPositionDefocusSA = [fixedX, fixedY, fixedDefocus, optimizedInclination, fixedAzimuth];
        % end



        function currentlnpdf = lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus) 
            currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
            currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - log(gamma(z+1)) , 'all');
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

