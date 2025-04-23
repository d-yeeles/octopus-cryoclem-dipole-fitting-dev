classdef FitPSF_ML_reparam2
    
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
                obj.estimatesPositionDefocus = fitting(obj);
            end
        end
        
        %% Fit
        function estimatesPositionDefocus = fitting(obj)
            parPsfEstimate = FitPSF_ML_reparam2.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(0, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift;

            % !!! don't hard-code these, inherit from elsewhere !!!
            parPsfEstimate.nPixels = 19;
            parPsfEstimate.pixelSize = Length(52, 'nm');
            parPsfEstimate.wavelength = Length(500, 'nm');
            parPsfEstimate.refractiveIndices(2) = 2.17;
            parPsfEstimate.refractiveIndices(1) = 1.31;
            parPsfEstimate.objectiveNA = 2.17;
            parPsfEstimate.refractiveIndices = [1.31 , 2.17 , 2.17];
            parPsfEstimate.pixelSensitivityMask = PixelSensitivity.uniform(9);
            parPsfEstimate.nDiscretizationBFP = 129;

            % psfEstimate = PSF(parPsfEstimate);
            if strcmpi(obj.model, 'hinterer')
                % parPsfEstimate.pixelSensitivityMask = PixelSensitivity.uniform(9);
                psfEstimate = PSF(parPsfEstimate);
            elseif strcmpi(obj.model, 'mortensen')
                % parPsfEstimate.pixelSensitivityMask = PixelSensitivity.uniform(1);
                psfEstimate = PSF_mortensen(parPsfEstimate);
            elseif strcmpi(obj.model, 'gaussian')
                % parPsfEstimate.pixelSensitivityMask = PixelSensitivity.uniform(1);
                psfEstimate = PSF_gaussian(parPsfEstimate);
            else
                error('Unknown model type: %s', obj.model);
            end

            % dave apr 2025 - did hinterer need this?
            psfImage = obj.image;% ./ norm(obj.image);
            
            % plotObjectiveFunction(obj, psfEstimate, psfImage);

            % estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate);%, estimatesPositionDefocus.LS);

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

        function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate)
            % Define the negative log-likelihood function
            lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
            costFunction = @(x) -lnpdf(image, x); % Negate likelihood for minimization
        
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
        
            % Initialize best solutions
            bestAttempt = [];
            bestCost = Inf;
            bestTheta = NaN;
            
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

    %     function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate)%, startValues)
    %         % Define the negative log-likelihood function
    %         lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
    %         costFunction = @(x) -lnpdf(image, x); % Negate likelihood for minimization
    % 
    %         % constain solutions so that cos(az)^2 + sin(az)^2 = 1
    %         % nonlcon = @(x) deal(x(4)^2 + x(5)^2 + x(6)^2 - 1, []);
    %         nonlcon = @(x) deal([], x(4)^2 + x(5)^2 + x(6)^2 - 1);
    % 
    %         startValues = [
    %             obj.parameterStartValues.x.inNanometer, ...
    %             obj.parameterStartValues.y.inNanometer, ...
    %             obj.parameterStartValues.defocus.inNanometer, ...
    %             obj.parameterStartValues.newangle1, ...
    %             obj.parameterStartValues.newangle2, ... 
    %             obj.parameterStartValues.newangle3, ... 
    %             obj.parameterStartValues.photons, ... 
    %             ];
    % 
    %         % Define parameter bounds
    %         xBounds = obj.parameterBounds.x.inNanometer;
    %         yBounds = obj.parameterBounds.y.inNanometer;
    %         defocusBounds = obj.parameterBounds.defocus.inNanometer;
    %         newangle1Bounds = obj.parameterBounds.newangle1;
    %         newangle2Bounds = obj.parameterBounds.newangle2;
    %         newangle3Bounds = obj.parameterBounds.newangle3;
    %         photonsBounds = obj.parameterBounds.photons;
    % 
    %         lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), newangle1Bounds(1), newangle2Bounds(1), newangle3Bounds(1), photonsBounds(1)];
    %         upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), newangle1Bounds(2), newangle2Bounds(2), newangle3Bounds(2), photonsBounds(2)];
    % 
    %         % options = optimoptions(@fmincon, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
    %         options = optimoptions(@fmincon, ...
    %             'Display', 'off');%, ...
    %             % 'StepTolerance', 1e-15, ...          % Stop when the step size is less than 1e-6
    %             % 'OptimalityTolerance', 1e-15, ...    % Stop when the gradient is less than 1e-6
    %             % 'FunctionTolerance', 1e-15, ...        % Stop if the function value change is less than 1e-6
    %             % 'MaxIterations', 1e15, ...          % Stop after 1000 iterations
    %             % 'MaxFunctionEvaluations', 1e15 ... % Stop after 5000 function evaluations
    %             % );
    % 
    %         % first optimisation attempt
    %         firstAttempt = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %         firstCost = costFunction(firstAttempt);
    %         firstThetaEst = acos(firstAttempt(6));
    %         firstThetaEst = mod(firstThetaEst, pi/2);
    % 
    %         % Store best result so far
    %         bestAttempt = firstAttempt;
    %         bestCost = firstCost;
    %         bestTheta = firstThetaEst;
    % 
    %         % Second attempt: random
    % 
    %         fprintf('   Attempt 2\n');
    %         newStartValues = startValues;
    %         for j = 1:7
    %             newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
    %         end
    % 
    %         % Run optimization
    %         secondAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %         secondCost = costFunction(secondAttempt);
    %         secondThetaEst = acos(secondAttempt(6));
    %         secondThetaEst = mod(secondThetaEst, pi/2);
    % 
    %         % Update best if this is better
    %         if secondCost < bestCost
    %             bestAttempt = secondAttempt;
    %             bestCost = secondCost;
    %             bestTheta = secondThetaEst;
    %         end
    % 
    %         % Third attempt: random
    % 
    %         fprintf('   Attempt 3\n');
    %         newStartValues = startValues;
    %         for j = 1:7
    %             newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
    %         end
    % 
    %         % Run optimization
    %         thirdAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %         thirdCost = costFunction(thirdAttempt);
    %         thirdThetaEst = acos(thirdAttempt(6));
    %         thirdThetaEst = mod(thirdThetaEst, pi/2);
    % 
    %         % Update best if this is better
    %         if thirdCost < bestCost
    %             bestAttempt = thirdAttempt;
    %             bestCost = thirdCost;
    %             bestTheta = thirdThetaEst;
    %         end
    % 
    %         % If it really thinks theta=90, then try shifting it
    % 
    %         condition = abs(bestTheta - pi/2) >= 1*(pi/180);
    % 
    %         if condition
    %             estimatesPositionDefocusML = bestAttempt;
    %         else
    % 
    %             % Fourth attempt: random, but shift theta by pi/2
    % 
    %             fprintf('   Attempt 4\n');
    %             newStartValues = startValues;
    %             for j = 1:7
    %                 newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
    %             end
    %             % but use the best theta so far minus pi/2
    %             newStartValues(6) = bestAttempt(6) + 1;
    % 
    %             % Run optimization
    %             fourthAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %             fourthCost = costFunction(fourthAttempt);
    %             fourthThetaEst = acos(fourthAttempt(6));
    %             fourthThetaEst = mod(fourthThetaEst, pi/2);
    % 
    %             % Update best if this is better
    %             if fourthCost < bestCost
    %                 bestAttempt = fourthAttempt;
    %                 bestCost = fourthCost;
    %                 bestTheta = fourthThetaEst;
    %             end
    % 
    %             % Fifth attempt: random, but shift theta by pi/4
    % 
    %             fprintf('   Attempt 5\n');
    %             newStartValues = startValues;
    %             for j = 1:7
    %                 newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
    %             end
    %             % but use the best theta so far minus pi/4
    %             newStartValues(6) = bestAttempt(6) + 1/sqrt(2);
    % 
    %             % Run optimization
    %             fifthAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %             fifthCost = costFunction(fifthAttempt);
    %             fifthThetaEst = acos(fifthAttempt(6));
    %             fifthThetaEst = mod(fifthThetaEst, pi/2);
    % 
    %             % Update best if this is better
    %             if fifthCost < bestCost
    %                 bestAttempt = fifthAttempt;
    %                 bestCost = fifthCost;
    %                 bestTheta = fifthThetaEst;
    %             end
    % 
    %             % If it really REALLY thinks theta=90, then do a bunch more
    %             % until it thinks it's not (or you hit max number of attempts)
    % 
    %             condition = abs(bestTheta - pi/2) >= 0.001*(pi/180);
    % 
    %             if condition
    %                 estimatesPositionDefocusML = bestAttempt;
    %             else
    %                 % Keep trying until condition is met or max_attempts is reached
    %                 max_attempts = 10;
    %                 current_attempt = 6; % We've already done 5 attempts
    % 
    %                 while ~condition && current_attempt <= max_attempts
    %                     fprintf('   Attempt %d\n', current_attempt);
    %                     newStartValues = startValues;
    %                     for j = 1:7
    %                         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
    %                     end
    % 
    %                     % Run optimization
    %                     nextAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %                     nextCost = costFunction(nextAttempt);
    %                     nextThetaEst = acos(nextAttempt(6));
    %                     nextThetaEst = mod(nextThetaEst, pi/2);
    % 
    %                     % Update best if this is better
    %                     if nextCost < bestCost
    %                         bestAttempt = nextAttempt;
    %                         bestCost = nextCost;
    %                         bestTheta = nextThetaEst;
    % 
    %                         % Check if condition is now met
    %                         condition = abs(bestTheta - pi/2) >= 0.001*(pi/180);
    %                     end
    % 
    %                     current_attempt = current_attempt + 1;
    %                 end
    %             end
    % 
    % 
    % 
    % 
    %             % Last resort, repeat
    %             % until it thinks it's not (or you hit max number of attempts)
    % 
    %             condition = bestTheta - pi/2 > 0;
    % 
    %             if condition
    %                 estimatesPositionDefocusML = bestAttempt;
    %             else
    %                 % Keep trying until condition is met or max_attempts is reached
    %                 max_attempts = 20;
    % 
    %                 while ~condition && current_attempt <= max_attempts
    %                     fprintf('   Attempt %d\n', current_attempt);
    %                     newStartValues = startValues;
    %                     for j = 1:7
    %                         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
    %                     end
    % 
    %                     % Run optimization
    %                     nextAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
    %                     nextCost = costFunction(nextAttempt);
    %                     nextThetaEst = acos(nextAttempt(6));
    %                     nextThetaEst = mod(nextThetaEst, pi/2);
    % 
    %                     % Update best if this is better
    %                     if nextCost < bestCost
    %                         bestAttempt = nextAttempt;
    %                         bestCost = nextCost;
    %                         bestTheta = nextThetaEst;
    % 
    %                         % Check if condition is now met
    %                         condition = abs(bestTheta - pi/2) >= 0.001*(pi/180);
    %                     end
    % 
    %                     current_attempt = current_attempt + 1;
    %                 end
    % 
    %             end
    % 
    %             estimatesPositionDefocusML = bestAttempt;
    %         end
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    %         end
    % 
    %         % add in something to save the value of the objective function
    %         estimatesPositionDefocusML(end+1) = bestCost;
    % 
    % end













            % % % without 180 fudge
            % % % estimatesPositionDefocusML = fminunc(@(x) -lnpdf(image, x), startValues, options);
            % % estimatesPositionDefocusML = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
            % % % add in something to save the value of the objective function
            % % estimatesPositionDefocusML(end+1) = costFunction(estimatesPositionDefocusML);
            % 
            % 
            % % with 180 fudge
            % % first optimisation attempt
            % firstAttempt = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
            % firstCost = costFunction(firstAttempt);
            % 
            % % Store best result so far
            % bestAttempt = firstAttempt;
            % bestCost = firstCost;
            % 
            % % Configuration for multiple starts
            % numAdditionalAttempts = 4;  % You were doing 4 more attempts
            % attempts = cell(1, numAdditionalAttempts + 1);
            % costs = zeros(1, numAdditionalAttempts + 1);
            % 
            % % Store first attempt
            % attempts{1} = firstAttempt;
            % costs(1) = firstCost;
            % 
            % % Check if we should try more optimizations
            % thetaEst = acos(firstAttempt(6));
            % phiEst = atan2(firstAttempt(5), firstAttempt(4));
            % thetaEst = mod(thetaEst, pi/2);
            % phiEst = mod(phiEst, 2*pi);
            % 
            % % Only try more if needed - use your condition
            % condition1 = abs(thetaEst - pi/2) >= 9999999999999; % your original condition
            % 
            % if condition1
            %     estimatesPositionDefocusML = firstAttempt; % this is good enough
            % else
            %     % disp('trying more')
            % 
            %     % Run additional optimization attempts with randomized starting points
            %     for i = 1:numAdditionalAttempts
            % 
            %         fprintf('   Attempt %d/%d\n', i+1, numAdditionalAttempts+1);
            % 
            %         % Create randomized starting values
            %         newStartValues = startValues;
            % 
            %         % % Randomize position and angles - consolidated from your various attempts
            %         % newStartValues(1) = (0 + 50 * (rand() - 0.5));  % x position
            %         % newStartValues(2) = (0 + 50 * (rand() - 0.5));  % y position
            %         % newStartValues(3) = 0;  % defocus
            % 
            %         for j = 1:7
            %             newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         end
            % 
            %         % if i == 1
            %         %     for j = 1:3
            %         %         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         %     end
            %         %     newStartValues(4) = 1;
            %         %     newStartValues(5) = 0;
            %         %     newStartValues(6) = lowerBounds(6) + (upperBounds(6) - lowerBounds(6)) * rand;
            %         % 
            %         % elseif i == 2
            %         %     for j = 1:3
            %         %         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         %     end
            %         %     newStartValues(4) = 0.5;
            %         %     newStartValues(5) = 0.866;
            %         %     newStartValues(6) = lowerBounds(6) + (upperBounds(6) - lowerBounds(6)) * rand;
            %         % 
            %         % elseif i == 3
            %         %     for j = 1:3
            %         %         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         %     end
            %         %     newStartValues(4) = -0.5;
            %         %     newStartValues(5) = 0.866;
            %         %     newStartValues(6) = lowerBounds(6) + (upperBounds(6) - lowerBounds(6)) * rand;
            %         % 
            %         % elseif i == 4
            %         %     for j = 1:3
            %         %         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         %     end
            %         %     newStartValues(4) = -1;
            %         %     newStartValues(5) = 0;
            %         %     newStartValues(6) = lowerBounds(6) + (upperBounds(6) - lowerBounds(6)) * rand;
            %         % 
            %         % elseif i == 5
            %         %     for j = 1:3
            %         %         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         %     end
            %         %     newStartValues(4) = -0.5;
            %         %     newStartValues(5) = -0.866;
            %         %     newStartValues(6) = lowerBounds(6) + (upperBounds(6) - lowerBounds(6)) * rand;
            %         % 
            %         % elseif i == 6
            %         %     for j = 1:3
            %         %         newStartValues(j) = lowerBounds(j) + (upperBounds(j) - lowerBounds(j)) * rand;
            %         %     end
            %         %     newStartValues(4) = 0.5;
            %         %     newStartValues(5) = -0.866;
            %         %     newStartValues(6) = lowerBounds(6) + (upperBounds(6) - lowerBounds(6)) * rand;
            %         % 
            %         % end
            % 
            %         % Run optimization
            %         currentAttempt = fmincon(@(x) -lnpdf(image, x), newStartValues, [], [], [], [], lowerBounds, upperBounds, nonlcon, options);
            %         currentCost = costFunction(currentAttempt);
            % 
            %         % Store result
            %         attempts{i+1} = currentAttempt;
            %         costs(i+1) = currentCost;
            % 
            %         % Update best if this is better
            %         if currentCost < bestCost
            %             bestAttempt = currentAttempt;
            %             bestCost = currentCost;
            %         end
            % 
            % 
            %     end
            % 
            %     % % Find the attempt with minimum cost
            %     % [~, minIndex] = min(costs);
            %     estimatesPositionDefocusML = bestAttempt;%attempts{minIndex};
            % end
            % 
            % % % Only doing x and y
            % % estimatesPositionDefocusML = firstAttempt;
            % % estimatesPositionDefocusML(end+1) = obj.parameterStartValues.defocus.inNanometer;
            % % estimatesPositionDefocusML(end+1) = obj.parameterStartValues.newangle1;
            % % estimatesPositionDefocusML(end+1) = obj.parameterStartValues.newangle2;
            % % estimatesPositionDefocusML(end+1) = obj.parameterStartValues.newangle3;
            % % estimatesPositionDefocusML(end+1) = obj.parameterStartValues.photons;
            % 
            % % add in something to save the value of the objective function
            % estimatesPositionDefocusML(end+1) = bestCost;%costFunction(estimatesPositionDefocusML);
        % 
        % end

        function currentlnpdf = lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus) 
            currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
            %currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - log(gamma(z+1)) , 'all');
            % dave feb 2025 - replacing with an all-in-one log gamma function just in case faster
            currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - gammaln(z+1) , 'all');
        end

        
        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)

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

            % currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
            % for k=1:size(psfEstimate.stageDrift.motion,1)
            %     aberrationCoeffs = getAberrations(psfEstimate,k);
            %     fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
            %     currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
            % end

            % ----------
            % dave jan 2025
            % doing more than the reduced form they were doing            
            % disp(xEstimate)
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

                % Show result
                % disp(currentPsf)

                currentPsf = double(currentPsf);

                % this bit is the equivalent of the getIntensitiesCamera()
                % stuff done in hinterer
                totalIntensity = sum(sum(currentPsf));
                currentPsf = currentPsf / totalIntensity * photonEstimate;
                currentFitPSF = currentPsf;

                % 
                % % !!! don't just hard-code this, read it in from elsewhere !!!
                % parTestPsf.nPixels = 19;
                % parTestPsf.wavelength = Length(500,'nm');
                % parTestPsf.objectiveNA = 2.17;
                % parTestPsf.objectiveFocalLength = Length(770,'mu');
                % parTestPsf.refractiveIndices = [1.31 2.17 2.17]; % [RI_specimen, RI_intermed, RI_immoil]
                % parTestPsf.nDiscretizationBFP = 129;%129; % change to like 501 or 1001 to make the ground truth dots line up better with the image
                % parTestPsf.pixelSize = Length(52,'nm');
                % parTestPsf.pixelSensitivityMask = PixelSensitivity.uniform(9);
                % parTestPsf.nPhotons = 2000;%1e10; % number of photons per spot - remember the image is made up by superimposing loads of individual images
                % parTestPsf.backgroundNoise = 0;
                % parTestPsf.shotNoise = 1;
                % parTestPsf.dipole = Dipole(inclination, azimuth); % dave jan 2025 - adding angle optimiser
                % parTestPsf.position = Length([xEstimate, yEstimate, 0], 'nm');
                % 
                % % parTestPsf = FitPSF_ML_reparam2.readParametersEstimate(obj.psf);
                % % parTestPsf.dipole = Dipole(inclination, azimuth); % dave jan 2025 - adding angle optimiser
                % % parTestPsf.position = Length([xEstimate, yEstimate, 0], 'nm');
                % % parTestPsf.nPhotons = obj.nPhotonEstimate;
                % % parTestPsf.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');
                % % parTestPsf.backgroundNoise = 0; % background noise is added later
                % % parTestPsf.pixelSensitivityMask = obj.pixelSensitivityMask;
                % % parTestPsf.stageDrift = obj.stageDrift;
                % % parTestPsf.nPixels = 19; % !!! don't hard-code this, inherit it from elsewhere !!!
                % test_psf = PSF_mortensen(parTestPsf);
                % 
                % % Output as tif
                % counter = 1000*rand();
                % output_path = sprintf('/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulate-multiple/output/1spot_allmodels/sims_mortensen_45/fitting/sim_frame%06d.tif', round(counter));
                % psf_total_image = uint32(currentPsf);
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

                % disp(['Raw Python PSF sum: ', num2str(sum(currentPsf,'all'))]);

                % currentPsf = ones(psfEstimate.nPixels,psfEstimate.nPixels); 

                % totalIntensity = sum(currentPsf,'all');
                % currentPsf = currentPsf ./ totalIntensity * photonEstimate + obj.noiseEstimate;
                % currentFitPSF = currentPsf ./ norm(currentPsf);

                % currentPsf = currentPsf + obj.noiseEstimate;
                % currentFitPSF = currentPsf;% ./ norm(currentPsf);
                % disp(['Processed PSF sum: ', num2str(sum(currentFitPSF,'all'))]);

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


                % ----------
    
                % photonEstimate = 1e9;
    
                % totalIntensity = sum(currentPsf,'all');
                % currentPsf = currentPsf ./ totalIntensity * photonEstimate + obj.noiseEstimate;
                % currentFitPSF = currentPsf ./ norm(currentPsf);
    
                % dave jan 2025
                % disp(['X: ', num2str(lateralPositionAndDefocus(1))]);
                % disp(['Inclination: ', num2str(lateralPositionAndDefocus(4))]);
    
                % % dave jan 2025 - print to check if PSF is updated every time
                % persistent iterationCounter;  % Keeps track of iterations
                % if isempty(iterationCounter)
                %     iterationCounter = 0;
                % end
                % iterationCounter = iterationCounter + 1;  % Increment counter
                % if mod(iterationCounter, 100) == 0  % Only output on every 10th iteration
                %     outputDirectory = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/optimiser_output';  % Define the output directory
                %     if ~exist(outputDirectory, 'dir')
                %         mkdir(outputDirectory);  % Create the directory if it doesn't exist
                %     end
                %     timestamp = datestr(now, 'yyyymmdd_HHMMSS_FFF'); 
                %     filename = sprintf('iteration_%s.png', timestamp);
                %     imwrite(mat2gray(currentFitPSF), fullfile(outputDirectory, filename));
                % end

        end

    end

    methods (Static)
        par = readParametersEstimate(psf);
    end
end