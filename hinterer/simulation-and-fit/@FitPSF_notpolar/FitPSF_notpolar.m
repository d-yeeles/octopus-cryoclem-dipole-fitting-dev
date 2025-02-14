classdef FitPSF_notpolar
    
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
        function obj = FitPSF_notpolar(psf, par)
            if nargin > 1
                obj = setInputParameters('FitPSF_notpolar', obj, par);
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
            parPsfEstimate = FitPSF_notpolar.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(obj.angleInclinationEstimate, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift; 

            psfEstimate = PSF(parPsfEstimate);

            psfImage = obj.image ./ norm(obj.image);
            
            plotObjectiveFunction(obj, psfEstimate, psfImage);%, 0);  
            % plotObjectiveFunction(obj, psfEstimate, psfImage, pi/4);  
            % plotObjectiveFunction(obj, psfEstimate, psfImage, pi/2);  
            % plotObjectiveFunction(obj, psfEstimate, psfImage, 3*pi/4);  
            % plotObjectiveFunction(obj, psfEstimate, psfImage, pi);  

            estimatesPositionDefocus.LS = fitLeastSquaresPSF(obj, psfImage, psfEstimate);
            estimatesPositionDefocus.ML = fitMaxLikelihoodPSF(obj, psfImage, psfEstimate, estimatesPositionDefocus.LS);

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

        function estimatesPositionDefocusML = fitMaxLikelihoodPSF(obj, image, psfEstimate, startValues)
            lnpdf = @(z,lateralPositionAndDefocus) lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus);
            options = optimoptions(@fmincon, 'Display', 'off', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
            % dave jan 2025
            % using constrained version, because why wouldnt you?
            defocusBounds = obj.parameterBounds.defocus.inNanometer;
            xBounds = obj.parameterBounds.x.inNanometer;
            yBounds = obj.parameterBounds.y.inNanometer;
            azimuthBounds = obj.parameterBounds.azimuth;
            lowerBounds = [xBounds(1), yBounds(1), defocusBounds(1), azimuthBounds(1)]; % unfixed azimuth
            upperBounds = [xBounds(2), yBounds(2), defocusBounds(2), azimuthBounds(2)]; % unfixed azimuth
            estimatesPositionDefocusML = fmincon(@(x) -lnpdf(image, x), startValues, [], [], [], [], lowerBounds, upperBounds, [], options);
        end

        function currentlnpdf = lnpdfFunction(obj,psfEstimate,z,lateralPositionAndDefocus) 
            currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
            currentlnpdf = sum(z.*log(currentPSF)  - currentPSF - log(gamma(z+1)) , 'all');
        end

        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)
            psfEstimate.position = Length([lateralPositionAndDefocus(1:2), 0], 'nm');
            psfEstimate.defocus = Length(lateralPositionAndDefocus(3), 'nm');
            lateralPositionAndDefocus(4) = mod(lateralPositionAndDefocus(4), 2*pi);

            psfEstimate.dipole = Dipole(obj.angleInclinationEstimate, lateralPositionAndDefocus(4)); % dave jan 2025 - adding angle optimiser

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

        end

        % % heat map objective function plot
        % function plotObjectiveFunction(obj, psfEstimate, image)
        %     % Define radial and azimuthal range
        %     rVals = linspace(0, 60, 60);
        %     azimuthVals = linspace(0, 2*pi, 180);
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
        %             % Convert (r, θ) to (x, y)
        %             x = rVals(i) * cos(azimuthVals(j));
        %             y = rVals(i) * sin(azimuthVals(j));
        %             params = [x, y, defocus, azimuthVals(j)];
        %             costMatrix(i, j) = -obj.lnpdfFunction(psfEstimate, image, params);
        %         end
        %     end
        % 
        %     % Convert polar to Cartesian for 3D surface plotting
        %     [AZI, R] = meshgrid(azimuthVals, rVals);
        %     [X, Y] = pol2cart(AZI, R);
        % 
        %     % Create a 3D surface plot
        %     figure;
        %     surf(X, Y, costMatrix, 'EdgeColor', 'none'); % Remove edges for smooth look
        %     view(2); % Top-down view
        %     colorbar;
        %     colormap jet;
        %     xlabel('X Position (nm)');
        %     ylabel('Y Position (nm)');
        %     title('Objective Function in Polar Coordinates (Defocus = 0)');
        % end

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


        % % 3D stack in x,y,az
        % function plotObjectiveFunction(obj, psfEstimate, image)
        %     % Define the grid of (x, y, azimuth) points
        %     numPoints = 20;  % Adjust for finer resolution
        %     xVals = linspace(-1000, 1000, 20);
        %     yVals = linspace(-1000, 1000, 20);
        %     azimuthVals = linspace(0, 2*pi, 20);
        % 
        %     % Initialize storage for results
        %     [X, Y, Azimuth] = meshgrid(xVals, yVals, azimuthVals);
        %     objFuncValues = zeros(size(X));
        % 
        %     % Evaluate the objective function at each sampled point
        %     for i = 1:numel(X)
        %         testPoint = [X(i), Y(i), 0, Azimuth(i)]; % Defocus fixed at 0
        %         objFuncValues(i) = -obj.lnpdfFunction(psfEstimate, image, testPoint); % Negative because we minimize
        %     end
        % 
        %     % 3D scatter plot
        %     figure;
        %     scatter3(X(:), Y(:), Azimuth(:), 50, objFuncValues(:), 'filled');
        %     colorbar;
        %     xlabel('X Position (nm)');
        %     ylabel('Y Position (nm)');
        %     zlabel('Azimuth (rad)');
        %     title('Objective Function Landscape');
        %     colormap jet;  % Use a colorful colormap
        %     view(135, 30); % Adjust viewing angle
        % end

        
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

