classdef FitPSF

    % dave december 2024
    % change startValues etc. to accept orientation estimates

    properties
        psf PSF 
        image
        noiseEstimate (1,1) {mustBeNonnegative, mustBeGreaterThanOrEqual(noiseEstimate, 1e-5)} = 1e-5  % >= 1e-5 for numeric stability of log(psf)
        nPhotonEstimate (1,1) {mustBeNonnegative} 
        stageDrift StageDrift = NoStageDrift()

        pixelSensitivityMask = PixelSensitivity.uniform(3)
        
        parameterBounds = struct( ...
            'x',Length([-800 800],'nm'), ...
            'y',Length([-800 800],'nm'), ...
            'defocus',Length([-2000 2000],'nm'), ...
            'inclination', [0 pi], ...
            'azimuth', [0 pi])

        parameterStartValues = struct( ...
            'x',Length(-100+200*rand(),'nm'), ...
            'y',Length(-100+200*rand(),'nm'), ...
            'defocus', Length(-500+1000*rand(),'nm'), ...
            'inclination', 0, ...
            'azimuth', 0)
        
        % Fit result
        estimatesPositionDefocus
    end
    
    methods
        function obj = FitPSF(psf, par)
            if nargin > 1
                obj = setInputParameters('FitPSF', obj, par);
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
            parPsfEstimate = FitPSF.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(0, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift; 

            psfEstimate = PSF(parPsfEstimate);

            psfImage = obj.image ./ norm(obj.image);
            
            estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);
            % estimatesPositionDefocus.GS = fitGlobalSearchPSF(obj, psfImage, psfEstimate);
        end
        
        function estimatesPositionDefocusLS = fitLeastSquaresPSF(obj, image, psfEstimate)
            funPsf = @(lateralPositionAndDefocus,xdata) createFitPSF(obj, psfEstimate, lateralPositionAndDefocus);
            xdata = zeros(obj.psf.nPixels,obj.psf.nPixels);
            options = optimoptions('lsqcurvefit','Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 5e-7, 'Display', 'off');  % Custom output function);
            
            startValues = [obj.parameterStartValues.x.inNanometer, ...
                obj.parameterStartValues.y.inNanometer, ...
                obj.parameterStartValues.defocus.inNanometer, ...
                obj.parameterStartValues.inclination, ...
                obj.parameterStartValues.azimuth];
            
            defocusBounds = obj.parameterBounds.defocus.inNanometer;
            xBounds = obj.parameterBounds.x.inNanometer;
            yBounds = obj.parameterBounds.y.inNanometer;
            inclinationBounds = obj.parameterBounds.inclination;
            azimuthBounds = obj.parameterBounds.azimuth;
            lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), inclinationBounds(1), azimuthBounds(1)];
            upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), inclinationBounds(2), azimuthBounds(2)];
            
            estimatesPositionDefocusLS = lsqcurvefit(funPsf, startValues, xdata, image, lowerBounds, upperBounds, options);

            disp(funPsf)
            disp(estimatesPositionDefocusLS)

        end
        
        function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate, startValues)
            lnpdf = @(z,lateralPositionAndDefocus) lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus);
            options = optimoptions(@fminunc, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
            estimatesPositionDefocusML = fminunc(@(x) -lnpdf(image,x), startValues, options);
            disp(estimatesPositionDefocusML)
        end

        % dave december 2024 - optimise over angles
        % function estimatesOrientationML = fitMaxLikelihoodPSF(obj, image, psfEstimate, startValues)
        %     lnpdf = @(z,lateralPositionAndDefocus) lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus);
        %     options = optimoptions(@fminunc, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
        %     estimatesOrientationML = fminunc(@(x) -lnpdf(image,x), startValues, options);
        % end

        % % dave december 2024 - this and all references to GS, GlobalSearch
        % function estimatesPositionDefocusGS = fitGlobalSearchPSF(obj, image, psfEstimate)
        % 
        %     % Define the log-likelihood function
        %     lnpdf = @(z, lateralPositionAndDefocus) lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus);
        % 
        %     % Define xdata (as in the original code, but this could vary)
        %     xdata = zeros(obj.psf.nPixels, obj.psf.nPixels);
        % 
        %     % Define start values from the object (same as in the original code)
        %     startValues = [obj.parameterStartValues.x.inNanometer, ...
        %         obj.parameterStartValues.y.inNanometer, ...
        %         obj.parameterStartValues.defocus.inNanometer];
        % 
        %     % Define the bounds for the parameters (from the object)
        %     defocusBounds = obj.parameterBounds.defocus.inNanometer;
        %     xBounds = obj.parameterBounds.x.inNanometer;
        %     yBounds = obj.parameterBounds.y.inNanometer;
        %     lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1)];
        %     upperBounds = [xBounds(2), yBounds(2), defocusBounds(2)];
        % 
        %     % Define the objective function (negative log-likelihood to minimize)
        %     objective = @(x) -lnpdf(image, x);  % Note the negative sign, because we want to minimize
        % 
        %     % Set up the optimization problem using fmincon
        %     problem = createOptimProblem('fmincon', ...
        %         'objective', objective, ...
        %         'x0', startValues, ...  % Initial guess
        %         'lb', lowerBounds, ...  % Lower bounds
        %         'ub', upperBounds);     % Upper bounds
        % 
        %     % Create a GlobalSearch object
        %     gs = GlobalSearch('NumTrialPoints', 500, 'Display', 'iter');  % Adjust trial points as needed
        % 
        %     % Run the optimization with GlobalSearch and fmincon
        %     estimatesPositionDefocusGS = run(gs, problem);
        % end


        function currentlnpdf = lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus) 
            currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
            currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - log(gamma(z+1)) , 'all');
        end

        
        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)
            psfEstimate.position = Length([lateralPositionAndDefocus(1:2), 0], 'nm');
            psfEstimate.defocus = Length(lateralPositionAndDefocus(3), 'nm');
            psfEstimate.dipole = Dipole(lateralPositionAndDefocus(4), lateralPositionAndDefocus(5));
            currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
            for k=1:size(psfEstimate.stageDrift.motion,1)
                aberrationCoeffs = getAberrations(psfEstimate,k);
                fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
                currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
            end
            totalIntensity = sum(currentPsf,'all');
            currentPsf = currentPsf ./ totalIntensity * obj.nPhotonEstimate + obj.noiseEstimate;
            currentFitPSF = currentPsf ./ norm(currentPsf);
        end

        
        %% Plot
        function plot(obj)
            plot(obj.psf)
            hold on
            size = 12;
            width = 2.5;
            center = (obj.psf.nPixels+1)/2;
            % plot(center+obj.estimatesPositionDefocus.LS(1)/100, center+obj.estimatesPositionDefocus.LS(2)/100,'Marker','o','MarkerSize',size,'Color','black','LineWidth', width)
            plot(center+obj.estimatesPositionDefocus.ML(1)/100, center+obj.estimatesPositionDefocus.ML(2)/100,'Marker','+','MarkerSize',size,'Color',[1 1 1]*0.8,'LineWidth', width)
            % % plot(center+obj.estimatesPositionDefocus.GS(1)/100, center+obj.estimatesPositionDefocus.GS(2)/100,'Marker','square','MarkerSize',size,'Color','r','LineWidth', width)
            plot(obj.psf.positionInPixelFromOrigin(1), obj.psf.positionInPixelFromOrigin(2),'Marker','x','MarkerSize',size,'Color','red','LineWidth', width)
            axis equal
            axis tight
            cb = colorbar;
            ylabel(cb,'Intensity','FontSize',15)
        end
    end

    methods (Static)
        par = readParametersEstimate(psf);
    end
end

