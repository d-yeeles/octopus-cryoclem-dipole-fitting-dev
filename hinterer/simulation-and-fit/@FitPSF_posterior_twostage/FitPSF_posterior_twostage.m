classdef FitPSF_posterior_twostage

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
            'theta', [0, pi/2], ...          % Inclination angle (0 to π/2)
            'phi', [0, 2*pi], ...            % Azimuth angle (0 to 2π)
            'photons', [1, 1e10]);           % Photon count

        parameterStartValues = struct( ...
            'x', Length(-100 + 200 * rand(), 'nm'), ...
            'y', Length(-100 + 200 * rand(), 'nm'), ...
            'defocus', Length(-500 + 1000 * rand(), 'nm'), ...
            'theta', pi/2 * rand(), ...      % Random inclination [0, π/2]
            'phi', 2*pi * rand(), ...        % Random azimuth [0, 2π]
            'photons', 1000);                % Initial photon estimate

        % Fit result - now includes nested sampling results
        estimatesPositionDefocus
        logEvidence  % Log evidence from nested sampling
        posteriorSamples  % Full posterior samples
        parameterUncertainties  % Parameter uncertainties
        credibleIntervals  % Credible intervals

        % BFP type indicator
        model = 'hinterer'  % Default to Hinterer BackFocalPlane

        % Nested sampling parameters
        nestedSamplingParams = struct( ...
            'Nlive', 60, ...
            'tolerance', 0.5, ...
            'Nmcmc', 0, ...  % 0 for multinest, >0 for MCMC
            'plotPosteriors', false, ...
            'estimationMethod', 'mean')  % 'mean', 'median', or 'map'
    end

    methods
        function obj = FitPSF_posterior_twostage(psf, par, model, autoFit)

            if nargin < 4
                autoFit = true;  % Default behavior for backward compatibility
            end

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
                obj = setInputParameters('FitPSF_posterior_twostage', obj, par);
            end

            if nargin > 0
                obj.psf = psf;
                obj.image = psf.image;
                obj.nPhotonEstimate = round(sum(sum(obj.image - obj.noiseEstimate)));
                % Only run fitting if autoFit is true
                if autoFit
                    [obj.estimatesPositionDefocus, obj.logEvidence, obj.posteriorSamples, ...
                     obj.parameterUncertainties, obj.credibleIntervals] = fitting(obj, model);
                end
            end
        end

        %% Fit using Nested Sampling
        function [estimatesPositionDefocus, logEvidence, posteriorSamples, parameterUncertainties, credibleIntervals] = fitting(obj, model)
            parPsfEstimate = FitPSF_posterior_twostage.readParametersEstimate(obj.psf);
            parPsfEstimate.dipole = Dipole(0, 0);
            parPsfEstimate.position = Length([0 0 0], 'nm');
            parPsfEstimate.nPhotons = obj.nPhotonEstimate;
            parPsfEstimate.defocus = Length(0, 'nm');
            parPsfEstimate.backgroundNoise = 0; % background noise is added later
            parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
            parPsfEstimate.stageDrift = obj.stageDrift;

            % Create PSF estimate based on model
            if strcmpi(obj.model, 'hinterer')
                psfEstimate = PSF(parPsfEstimate);
            elseif strcmpi(obj.model, 'mortensen')
                psfEstimate = PSF_mortensen(parPsfEstimate);
            elseif strcmpi(obj.model, 'gaussian')
                psfEstimate = PSF_gaussian(parPsfEstimate);
            else
                error('Unknown model type: %s', obj.model);
            end

            % Prepare image
            if strcmpi(model, 'gaussian')
                psfImage = obj.image ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
            elseif strcmpi(model, 'hinterer')
                psfImage = obj.image;% ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
            elseif strcmpi(model, 'mortensen')
                psfImage = obj.image;% ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
            end

            % Run nested sampling
            [logEvidence, nest_samples, post_samples] = fitNestedSamplingPSF(obj, psfImage, psfEstimate, model);

            % Extract point estimates
            estimatesPositionDefocus.ML = FitPSF_posterior_twostage.extractPointEstimates(post_samples, obj.nestedSamplingParams.estimationMethod);

            % Store full results
            posteriorSamples = post_samples;

            % Get parameter uncertainties
            [parameterUncertainties, credibleIntervals] = FitPSF_posterior_twostage.getParameterUncertainties(post_samples, 0.95);
        end

        function [logZ, nest_samples, post_samples] = fitNestedSamplingPSF(obj, image, psfEstimate, model)
            fprintf('Starting nested sampling PSF fit...\n');

            % CRITICAL FIX: Ensure ALL parameters are scalar doubles
            Nlive = double(obj.nestedSamplingParams.Nlive);
            if ~isscalar(Nlive)
                Nlive = Nlive(1);  % Take first element if it's an array
            end

            tolerance = double(obj.nestedSamplingParams.tolerance);
            if ~isscalar(tolerance)
                tolerance = tolerance(1);
            end

            Nmcmc = double(obj.nestedSamplingParams.Nmcmc);
            if ~isscalar(Nmcmc)
                Nmcmc = Nmcmc(1);
            end

            % Additional validation to ensure they're finite scalars
            if ~isfinite(Nlive) || Nlive < 1
                Nlive = 50;  % Safe default
            end
            if ~isfinite(tolerance) || tolerance <= 0
                tolerance = 0.1;  % Safe default  
            end
            if ~isfinite(Nmcmc) || Nmcmc < 0
                Nmcmc = 0;  % Safe default
            end

            % Force to integers where needed
            Nlive = round(Nlive);
            Nmcmc = round(Nmcmc);

            fprintf('Validated parameters: Nlive=%d, tolerance=%.3f, Nmcmc=%d\n', Nlive, tolerance, Nmcmc);

            % % Ask user for MCMC iterations if not set
            % if isempty(Nmcmc) || Nmcmc == 0
            %     Nmcmc = input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
            %     if isempty(Nmcmc)
            %         Nmcmc = 0;
            %     end
            %     Nmcmc = round(double(Nmcmc));  % Ensure scalar integer
            % end
            Nmcmc = 0;
            Nmcmc = round(double(Nmcmc));

            % Prepare data for nested sampling - ensure it's properly formatted
            data = double(image);

            fprintf('Data validation:\n');
            fprintf('Data size: %s\n', mat2str(size(data)));
            fprintf('Data class: %s\n', class(data));
            fprintf('Data range: [%.2e, %.2e]\n', min(data(:)), max(data(:)));
            fprintf('Data finite: %s\n', mat2str(all(isfinite(data(:)))));

            % Define the likelihood function with additional validation
            likelihood = @(data_in, model_func, parnames, parvals) FitPSF_posterior_twostage.validateAndCallLikelihood(obj, psfEstimate, data_in, parnames, parvals);

            % Define the PSF model function (placeholder)
            model_func = @(data_in, parnames, parvals) ones(size(data_in));

            if strcmpi(model, 'gaussian')
                % Define priors for Gaussian model (4 parameters)
                xBounds = obj.parameterBounds.x.inNanometer;
                yBounds = obj.parameterBounds.y.inNanometer;
                defocusBounds = obj.parameterBounds.defocus.inNanometer;
                photonsBounds = obj.parameterBounds.photons;

                % Ensure all bounds are scalar doubles and finite
                xLower = double(xBounds(1));
                xUpper = double(xBounds(2));
                yLower = double(yBounds(1));
                yUpper = double(yBounds(2));
                defocusLower = double(defocusBounds(1));
                defocusUpper = double(defocusBounds(2));
                photonsLower = double(photonsBounds(1));
                photonsUpper = double(photonsBounds(2));

                % Additional safety checks
                if ~isfinite(xLower) || ~isfinite(xUpper) || xLower >= xUpper
                    xLower = -100; xUpper = 100;
                end
                if ~isfinite(yLower) || ~isfinite(yUpper) || yLower >= yUpper
                    yLower = -100; yUpper = 100;
                end
                if ~isfinite(defocusLower) || ~isfinite(defocusUpper) || defocusLower >= defocusUpper
                    defocusLower = -10; defocusUpper = 10;
                end
                if ~isfinite(photonsLower) || ~isfinite(photonsUpper) || photonsLower >= photonsUpper
                    photonsLower = 100; photonsUpper = 10000;
                end

                prior = {'x', 'uniform', xLower, xUpper, 'fixed'; ...
                         'y', 'uniform', yLower, yUpper, 'fixed'; ...
                         'defocus', 'uniform', defocusLower, defocusUpper, 'fixed'; ...
                         'photons', 'uniform', photonsLower, photonsUpper, 'fixed'};

            else
                % Define priors for Hinterer/Mortensen model (6 parameters: x, y, defocus, theta, phi, photons)
                xBounds = obj.parameterBounds.x.inNanometer;
                yBounds = obj.parameterBounds.y.inNanometer;
                defocusBounds = obj.parameterBounds.defocus.inNanometer;
                thetaBounds = obj.parameterBounds.theta;      % [0, π/2]
                phiBounds = obj.parameterBounds.phi;          % [0, 2π]
                photonsBounds = obj.parameterBounds.photons;

                % Ensure all bounds are scalar doubles and finite
                xLower = double(xBounds(1));
                xUpper = double(xBounds(2));
                yLower = double(yBounds(1));
                yUpper = double(yBounds(2));
                defocusLower = double(defocusBounds(1));
                defocusUpper = double(defocusBounds(2));
                thetaLower = double(thetaBounds(1));
                thetaUpper = double(thetaBounds(2));
                phiLower = double(phiBounds(1));
                phiUpper = double(phiBounds(2));
                photonsLower = double(photonsBounds(1));
                photonsUpper = double(photonsBounds(2));

                % Additional safety checking for all bounds
                bounds_check = [xLower, xUpper; yLower, yUpper; defocusLower, defocusUpper; ...
                               thetaLower, thetaUpper; phiLower, phiUpper; photonsLower, photonsUpper];

                for i = 1:size(bounds_check, 1)
                    if ~isfinite(bounds_check(i,1)) || ~isfinite(bounds_check(i,2)) || bounds_check(i,1) >= bounds_check(i,2)
                        switch i
                            case 1  % x bounds
                                xLower = -100; xUpper = 100;
                            case 2  % y bounds
                                yLower = -100; yUpper = 100;
                            case 3  % defocus bounds
                                defocusLower = -10; defocusUpper = 10;
                            case 4  % theta bounds
                                thetaLower = 0; thetaUpper = pi/2;
                            case 5  % phi bounds
                                phiLower = 0; phiUpper = 2*pi;
                            case 6  % photons bounds
                                photonsLower = 100; photonsUpper = 10000;
                        end
                    end
                end

                prior = {'x', 'uniform', xLower, xUpper, 'fixed'; ...
                         'y', 'uniform', yLower, yUpper, 'fixed'; ...
                         'defocus', 'uniform', defocusLower, defocusUpper, 'fixed'; ...
                         'theta', 'uniform', thetaLower, thetaUpper, 'fixed'; ...
                         'phi', 'uniform', phiLower, phiUpper, 'fixed'; ...
                         'photons', 'uniform', photonsLower, photonsUpper, 'fixed'};
            end

            % Final validation of prior structure - ensure all entries are scalar
            fprintf('Prior validation:\n');
            for i = 1:size(prior, 1)
                param_name = prior{i,1};
                prior_type = prior{i,2};
                lower_bound = prior{i,3};
                upper_bound = prior{i,4};

                % Force to scalar if needed
                if ~isscalar(lower_bound)
                    lower_bound = lower_bound(1);
                    prior{i,3} = lower_bound;
                end
                if ~isscalar(upper_bound)
                    upper_bound = upper_bound(1);
                    prior{i,4} = upper_bound;
                end

                fprintf('  %s: %s [%.6f, %.6f]\n', param_name, prior_type, lower_bound, upper_bound);

                if ~isscalar(lower_bound) || ~isscalar(upper_bound)
                    error('Prior bounds must be scalar for parameter %s', param_name);
                end
                if ~isfinite(lower_bound) || ~isfinite(upper_bound)
                    error('Prior bounds must be finite for parameter %s', param_name);
                end
                if lower_bound >= upper_bound
                    error('Lower bound must be less than upper bound for parameter %s', param_name);
                end
            end

            extraparams = {};  % Ensure this is empty cell array

            fprintf('\nRunning nested sampler with %d live points...\n', Nlive);

            % Run nested sampling with all scalar parameters
            try
                % CRITICAL: Ensure all optional parameters are also scalar
                diffevfrac_val = double(2);
                covfrac_val = double(10);
                walkfrac_val = double(0);
                stretchfrac_val = double(0);

                [logZ, nest_samples, post_samples] = nested_sampler(...
                    data, ...                    % data matrix
                    Nlive, ...                   % scalar integer  
                    tolerance, ...               % scalar double
                    likelihood, ...              % function handle
                    model_func, ...              % function handle
                    prior, ...                   % validated cell array
                    extraparams, ...             % empty cell array
                    'Nmcmc', Nmcmc, ...         % scalar integer
                    'diffevfrac', diffevfrac_val, ...  % scalar double
                    'covfrac', covfrac_val, ...        % scalar double
                    'walkfrac', walkfrac_val, ...      % scalar double
                    'stretchfrac', stretchfrac_val);   % scalar double

                fprintf('Nested sampling completed successfully. Log evidence: %.2f\n', logZ);

            catch ME
                fprintf('\n=== NESTED SAMPLER ERROR ===\n');
                fprintf('Error message: %s\n', ME.message);
                fprintf('Error identifier: %s\n', ME.identifier);

                for i = 1:length(ME.stack)
                    fprintf('Stack frame %d: %s (line %d)\n', i, ME.stack(i).file, ME.stack(i).line);
                end

                fprintf('============================\n');
                rethrow(ME);
            end

            % Optional: Plot posterior distributions  
            if obj.nestedSamplingParams.plotPosteriors
                fprintf('Plotting posterior distributions...\n');

                try
                    if strcmpi(model, 'gaussian')
                        for i = 1:4
                            posteriors(post_samples, i, {prior{i,1}});
                        end
                        posteriors(post_samples, [1 2], {prior{1,1}, prior{2,1}});
                    else
                        % For spherical coordinates (6 parameters)
                        for i = 1:6
                            posteriors(post_samples, i, {prior{i,1}});
                        end
                        posteriors(post_samples, [1 2], {prior{1,1}, prior{2,1}});
                        % Plot angles together
                        posteriors(post_samples, [4 5], {prior{4,1}, prior{5,1}});
                    end
                catch plot_ME
                    fprintf('Warning: Could not create posterior plots: %s\n', plot_ME.message);
                end
            end
        end

        function currentlnpdf = lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus) 
            try
                currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 

                % Ensure currentPSF and z have the same size
                if ~isequal(size(currentPSF), size(z))
                    error('Size mismatch: currentPSF size = %s, z size = %s', ...
                        mat2str(size(currentPSF)), mat2str(size(z)));
                end

                % Compute log-likelihood element-wise, then sum
                log_psf = log(currentPSF);

                % Check for invalid values
                if any(~isfinite(log_psf(:))) || any(~isfinite(currentPSF(:)))
                    currentlnpdf = -inf;
                    return;
                end

                % Compute Poisson log-likelihood
                individual_terms = z .* log_psf - currentPSF - gammaln(z + 1);

                % Sum all terms to get scalar result
                currentlnpdf = sum(individual_terms(:));

                % Final check that result is scalar
                if ~isscalar(currentlnpdf)
                    error('lnpdfFunction returned non-scalar result');
                end

            catch ME
                fprintf('Error in lnpdfFunction: %s\n', ME.message);
                currentlnpdf = -inf;
            end
        end

        function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)

            if strcmpi(obj.model, 'gaussian')
                % Gaussian model: x, y, defocus, photons
                xEstimate = lateralPositionAndDefocus(1);
                yEstimate = lateralPositionAndDefocus(2);
                photonEstimate = lateralPositionAndDefocus(4);

                psfEstimate.position = Length([xEstimate, yEstimate, 0], 'nm');
                psfEstimate.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');

            else
                % Non-Gaussian models: x, y, defocus, theta, phi, photons
                xEstimate = lateralPositionAndDefocus(1);
                yEstimate = lateralPositionAndDefocus(2);
                % defocusEstimate = lateralPositionAndDefocus(3);  % Currently not used
                theta = lateralPositionAndDefocus(4);   % Inclination angle
                phi = lateralPositionAndDefocus(5);     % Azimuth angle  
                photonEstimate = lateralPositionAndDefocus(6);

                % theta and phi are already in the correct format
                inclination = theta;
                azimuth = phi;

                % Ensure angles are in valid ranges
                inclination = mod(inclination, pi/2);    % [0, π/2]
                azimuth = mod(azimuth, 2*pi);            % [0, 2π]

                psfEstimate.position = Length([xEstimate, yEstimate, 0], 'nm');
                psfEstimate.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');
                psfEstimate.dipole = Dipole(inclination, azimuth);
            end

            % dave jan 2025
            % doing more than the reduced form they were doing            
            if strcmpi(obj.model, 'mortensen')

                % Use Mortensen model

                % Call Python stuff
                % Try two specific locations for the Python module (local vs cluster)
                possibleDirs = {
                    fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
                    '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits'
                };

                % Try each location until we find the Python file
                pyModuleFound = false;
                for i = 1:length(possibleDirs)
                    pyDir = possibleDirs{i};

                    % Check if directory exists
                    if ~exist(pyDir, 'dir')
                        continue;  % Skip to next directory
                    end

                    % Check if Python file exists in this directory
                    pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
                    if exist(pyFilePath, 'file')

                        % Add to Python path
                        if count(py.sys.path(), pyDir) == 0
                            py.sys.path().insert(int32(0), pyDir);
                        end

                        pyModuleFound = true;
                        break;  % Found the file, stop looking
                    end
                end

                % If the module wasn't found in any location, show an error
                if ~pyModuleFound
                    error(['Python module "vectorized_mortensen_flipy_phi_theta.py" not found in any of the specified locations. ', ...
                           'Please ensure the file exists in one of these directories:\n', ...
                           '- %s\n', ...
                           '- %s'], possibleDirs{1}, possibleDirs{2});
                end

                % Run the function defined in your python file
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

                % Convert to matlab array
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

                currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
                for k=1:size(psfEstimate.stageDrift.motion,1)
                    % Apply aberrations
                    aberrationCoeffs = getAberrations(psfEstimate,k);
                    fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
                    % Get image from BFP field
                    currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
                end

                totalIntensity = sum(currentPsf,'all');
                currentPsf = currentPsf ./ totalIntensity * photonEstimate + obj.noiseEstimate;
                currentFitPSF = currentPsf;% !!! maybe swap this back - 8th May 2025

            end

            currentFitPSF = adjustExcitation(psfEstimate, currentFitPSF);
            currentFitPSF = applyShotNoise(psfEstimate, currentFitPSF);
            currentFitPSF = addBackgroundNoise(psfEstimate, currentFitPSF);
        end

        % Method to set nested sampling parameters
        function obj = setNestedSamplingParams(obj, varargin)
            % Set nested sampling parameters
            % Usage: obj = obj.setNestedSamplingParams('Nlive', 1000, 'tolerance', 0.01, ...)

            for i = 1:2:length(varargin)
                param = varargin{i};
                value = varargin{i+1};

                if isfield(obj.nestedSamplingParams, param)
                    obj.nestedSamplingParams.(param) = value;
                else
                    warning('Unknown nested sampling parameter: %s', param);
                end
            end
        end

        function obj = runFittingWithParams(obj, varargin)
            % Set parameters first
            obj = obj.setNestedSamplingParams(varargin{:});

            % Then run fitting
            [obj.estimatesPositionDefocus, obj.logEvidence, obj.posteriorSamples, ...
             obj.parameterUncertainties, obj.credibleIntervals] = fitting(obj, obj.model);
        end

        % Method to get parameter summary
        function summary = getParameterSummary(obj)
            % Get a summary of parameter estimates and uncertainties

            if isempty(obj.posteriorSamples)
                error('No posterior samples available. Run fitting first.');
            end

            samples = obj.posteriorSamples.samples;

            if strcmpi(obj.model, 'gaussian')
                param_names = {'x', 'y', 'defocus', 'photons'};
            else
                param_names = {'x', 'y', 'defocus', 'theta', 'phi', 'photons'};
            end

            summary = struct();
            summary.parameter_names = param_names;
            summary.estimates = obj.estimatesPositionDefocus.ML(1:end-1);  % Exclude cost
            summary.uncertainties = obj.parameterUncertainties;
            summary.credible_intervals = obj.credibleIntervals;
            summary.log_evidence = obj.logEvidence;
            summary.estimation_method = obj.nestedSamplingParams.estimationMethod;

            % Display summary
            fprintf('\n=== Parameter Estimation Summary ===\n');
            fprintf('Model: %s\n', obj.model);
            fprintf('Estimation method: %s\n', summary.estimation_method);
            fprintf('Log evidence: %.2f\n', summary.log_evidence);
            fprintf('\nParameter estimates (±1σ, [95%% CI]):\n');

            for i = 1:length(param_names)
                fprintf('%s: %.2f ± %.2f [%.2f, %.2f]\n', ...
                    param_names{i}, ...
                    summary.estimates(i), ...
                    summary.uncertainties(i), ...
                    summary.credible_intervals(i,1), ...
                    summary.credible_intervals(i,2));
            end
            fprintf('\n');
        end
    end

    methods (Static)
        % function par = readParametersEstimate(psf)
        %     % dave apr 2025 - add this so we can inherit the right params when
        %     % fitting
        %     par = struct();
        %     par.nPixels = psf.nPixels;
        %     par.pixelSize = psf.pixelSize;
        %     par.wavelength = psf.wavelength;
        %     par.objectiveNA = psf.objectiveNA;
        %     par.pixelSensitivityMask = psf.pixelSensitivityMask;
        %     par.nDiscretizationBFP = psf.nDiscretizationBFP;
        % 
        %     % Fluorophore
        %     par.shotNoise = 0; 
        %     par.reducedExcitation = 0;
        % 
        %     % Microscope setup
        %     par.wavelength = psf.wavelength;
        %     if isprop(psf, 'astigmatism') || isfield(psf, 'astigmatism')
        %         par.astigmatism = psf.astigmatism;
        %     end
        %     par.objectiveNA = psf.objectiveNA;
        %     par.objectiveFocalLength = psf.objectiveFocalLength;
        %     par.refractiveIndices = psf.refractiveIndices;
        %     if isprop(psf, 'heightIntermediateLayer') || isfield(psf, 'heightIntermediateLayer')
        %         par.heightIntermediateLayer = psf.heightIntermediateLayer;
        %     end
        % 
        %     % Back focal plane
        %     if isprop(psf, 'phaseMask') || isfield(psf, 'phaseMask')
        %         par.phaseMask = psf.phaseMask;
        %     end
        %     par.nDiscretizationBFP = psf.nDiscretizationBFP;
        % 
        %     % Camera
        %     par.pixelSize = psf.pixelSize;
        % 
        %     % dave apr 2025
        %     par.nPixels = psf.nPixels;
        % end

        par = readParametersEstimate(psf);

        function logL = validateAndCallLikelihood(obj, psfEstimate, data, parnames, parvals)
            try
                % Simple validation
                if length(parnames) ~= length(parvals)
                    logL = -inf;
                    return;
                end

                % Call the nested likelihood function - note the argument order
                logL = FitPSF_posterior_twostage.lnpdfFunction_nested(obj, obj, psfEstimate, data, parnames, parvals);

                % Final validation
                if ~isscalar(logL) || ~isreal(logL)
                    fprintf('Warning: Likelihood function returned invalid result\n');
                    logL = -inf;
                end

            catch ME
                fprintf('Error in validateAndCallLikelihood: %s\n', ME.message);
                if length(ME.stack) > 0
                    fprintf('  at line %d in %s\n', ME.stack(1).line, ME.stack(1).file);
                end
                logL = -inf;
            end
        end

        function logL = lnpdfFunction_nested(obj, obj_ref, psfEstimate, data, parnames, parvals)
            % Convert parameter names and values to the format expected by your original function

            % Debug: Print what we're receiving
            persistent debug_count;
            if isempty(debug_count)
                debug_count = 0;
            end
            debug_count = debug_count + 1;

            try
                % Handle the case where nested_sampler might pass arguments differently
                % Sometimes nested_sampler doesn't follow the expected format exactly

                % If data is not a matrix, it might be that the arguments are in a different order
                if ~ismatrix(data) || ~isnumeric(data)
                    % Try to find the actual data in the arguments
                    % Sometimes the first argument might be something else
                    fprintf('Warning: data is not a numeric matrix, trying to find actual data...\n');

                    % Check if obj_ref contains the data we need
                    if isfield(obj_ref, 'image') && isnumeric(obj_ref.image)
                        actual_data = obj_ref.image;
                        fprintf('Found data in obj_ref.image: %s\n', mat2str(size(actual_data)));
                    elseif isnumeric(obj_ref) && ismatrix(obj_ref)
                        actual_data = obj_ref;
                        fprintf('Using obj_ref as data: %s\n', mat2str(size(actual_data)));
                    else
                        % Last resort: use the original data from obj
                        actual_data = obj.image;
                        fprintf('Using obj.image as fallback: %s\n', mat2str(size(actual_data)));
                    end
                else
                    actual_data = data;
                end

                % Ensure we have valid data
                if ~ismatrix(actual_data) || ~isnumeric(actual_data)
                    logL = -inf;
                    return;
                end

                % Build parameter vector in the expected order
                if length(parnames) == 4  % Gaussian model
                    params = zeros(1, 4);
                    for i = 1:length(parnames)
                        val = parvals{i};
                        if ~isscalar(val) || ~isreal(val) || ~isfinite(val)
                            logL = -inf;
                            return;
                        end

                        switch parnames{i}
                            case 'x'
                                params(1) = double(val);
                            case 'y'
                                params(2) = double(val);
                            case 'defocus'
                                params(3) = double(val);
                            case 'photons'
                                params(4) = double(val);
                        end
                    end
                else  % Full model with spherical angles (6 parameters)
                    params = zeros(1, 6);
                    for i = 1:length(parnames)
                        val = parvals{i};
                        if ~isscalar(val) || ~isreal(val) || ~isfinite(val)
                            logL = -inf;
                            return;
                        end

                        switch parnames{i}
                            case 'x'
                                params(1) = double(val);
                            case 'y'
                                params(2) = double(val);
                            case 'defocus'
                                params(3) = double(val);
                            case 'theta'
                                params(4) = double(val);
                            case 'phi'
                                params(5) = double(val);
                            case 'photons'
                                params(6) = double(val);
                        end
                    end
                end

                % Additional parameter validation
                if length(params) == 4 && params(4) <= 0  % Gaussian model photons
                    logL = -inf;
                    return;
                elseif length(params) == 6 && params(6) <= 0  % Full model photons
                    logL = -inf;
                    return;
                end

                % Call your original likelihood function with the correct object reference
                % Use obj (the FitPSF_posterior_twostage object) instead of obj_ref
                logL_result = obj.lnpdfFunction(psfEstimate, actual_data, params);

                % CRITICAL: Ensure the result is a scalar
                if ~isscalar(logL_result)
                    error('Likelihood function returned non-scalar result: size = %s', mat2str(size(logL_result)));
                end

                % Ensure log-likelihood is finite and real
                if ~isreal(logL_result) || ~isfinite(logL_result)
                    logL = -inf;
                else
                    logL = double(logL_result);  % Ensure it's a double scalar
                end

            catch ME
                % If any error occurs in likelihood computation, return very low likelihood
                fprintf('Error in likelihood computation: %s\n', ME.message);
                if length(ME.stack) > 0
                    fprintf('Error in file %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
                end
                logL = -inf;
            end
        end

        function model_output = psfModel_nested(obj, obj_ref, psfEstimate, data, parnames, parvals)
            % This function isn't strictly necessary for your PSF fitting since you're
            % computing likelihood directly from the image, but nested_sampler requires it
            model_output = ones(size(data));  % Placeholder
        end

        function estimatesPositionDefocusML = extractPointEstimates(post_samples, method)
            % Extract point estimates from posterior samples
            % 
            % Inputs:
            %   post_samples - posterior samples from nested sampling
            %   method - 'mean', 'median', 'map' (maximum a posteriori)
            %
            % Output:
            %   estimatesPositionDefocusML - parameter estimates in same format as original

            if nargin < 2
                method = 'mean';  % Default to posterior mean
            end

            % Handle different formats of post_samples
            if isstruct(post_samples) && isfield(post_samples, 'samples') && isfield(post_samples, 'logL')
                samples = post_samples.samples;
                logL_values = post_samples.logL;
            elseif iscell(post_samples) && length(post_samples) >= 2
                % Alternative format from nested_sampler
                samples = post_samples{1};
                logL_values = post_samples{2};
            else
                % Try to adapt to whatever format we have
                if isnumeric(post_samples)
                    % Assume post_samples is just the samples matrix and we don't have logL values
                    samples = post_samples;
                    logL_values = zeros(size(samples, 1), 1);  % Dummy values
                    fprintf('Warning: No log-likelihood values available. Using dummy values.\n');
                else
                    error('Unrecognized format for post_samples');
                end
            end

            switch lower(method)
                case 'mean'
                    % Posterior mean (most common choice)
                    estimatesPositionDefocusML = mean(samples, 1);

                case 'median'
                    % Posterior median (more robust to outliers)
                    estimatesPositionDefocusML = median(samples, 1);

                case 'map'
                    % Maximum a posteriori (MAP) estimate - sample with highest likelihood
                    [~, map_idx] = max(logL_values);
                    estimatesPositionDefocusML = samples(map_idx, :);

                otherwise
                    error('Unknown method. Use ''mean'', ''median'', or ''map''');
            end

            % Add the log-likelihood value at the end (equivalent to cost in original)
            % For MAP estimate, use the actual max likelihood
            % For mean/median, evaluate likelihood at that point
            if strcmpi(method, 'map') && ~isempty(logL_values)
                estimatesPositionDefocusML(end+1) = -logL_values(map_idx);  % Convert to cost (negative log-likelihood)
            else
                % For mean/median, we use the mean of the log-likelihood values
                if ~isempty(logL_values)
                    estimatesPositionDefocusML(end+1) = -mean(logL_values);
                else
                    estimatesPositionDefocusML(end+1) = 0;  % Default value if no log-likelihood available
                end
            end
        end

        function [param_uncertainties, credible_intervals] = getParameterUncertainties(post_samples, confidence_level)
            % Extract parameter uncertainties from posterior samples
            %
            % Inputs:
            %   post_samples - posterior samples from nested sampling  
            %   confidence_level - confidence level for credible intervals (default: 0.95)
            %
            % Outputs:
            %   param_uncertainties - standard deviations of parameters
            %   credible_intervals - credible intervals [lower, upper] for each parameter

            if nargin < 2
                confidence_level = 0.95;
            end

            % Handle different formats of post_samples
            if isstruct(post_samples) && isfield(post_samples, 'samples')
                samples = post_samples.samples;
            elseif iscell(post_samples) && length(post_samples) >= 1
                % Alternative format from nested_sampler
                samples = post_samples{1};
            else
                % Try to adapt to whatever format we have
                if isnumeric(post_samples)
                    % Assume post_samples is just the samples matrix
                    samples = post_samples;
                else
                    error('Unrecognized format for post_samples');
                end
            end

            % Parameter uncertainties (standard deviations)
            param_uncertainties = std(samples, 1);

            % Credible intervals
            alpha = 1 - confidence_level;
            lower_percentile = 100 * alpha / 2;
            upper_percentile = 100 * (1 - alpha / 2);

            credible_intervals = [prctile(samples, lower_percentile, 1); ...
                                 prctile(samples, upper_percentile, 1)]';
        end
    end
end



% %% unit vector version
% 
% classdef FitPSF_posterior_twostage
% 
%     properties
%         psf % PSF  % dave apr 2025 - removing this to avoid type coercion
%         image
%         % angleInclinationEstimate (1,1) {mustBeInFullRadialRange} % dave jan 2025 - commented out for adding angle optimiser
%         % angleAzimuthEstimate (1,1) {mustBeInFullRadialRange} % dave jan 2025 - commented out for adding angle optimiser
%         noiseEstimate (1,1) {mustBeNonnegative, mustBeGreaterThanOrEqual(noiseEstimate, 1e-5)} = 1e-5  % >= 1e-5 for numeric stability of log(psf)
%         nPhotonEstimate (1,1) {mustBeNonnegative} 
%         stageDrift StageDrift = NoStageDrift()
% 
%         pixelSensitivityMask = PixelSensitivity.uniform(3)
% 
%         parameterBounds = struct( ...
%             'x', Length([-800 800], 'nm'), ...
%             'y', Length([-800 800], 'nm'), ...
%             'defocus', Length([-2000 2000], 'nm'), ...
%             'unitvector1', [-1, 1], ...          % Inclination angle (0 to π/2)
%             'unitvector2', [-1, 1], ...            % Azimuth angle (0 to 2π)
%             'unitvector3', [0, 1], ...            % Azimuth angle (0 to 2π)
%             'photons', [1, 1e10]);           % Photon count
% 
%         parameterStartValues = struct( ...
%             'x', Length(-100 + 200 * rand(), 'nm'), ...
%             'y', Length(-100 + 200 * rand(), 'nm'), ...
%             'defocus', Length(-500 + 1000 * rand(), 'nm'), ...
%             'unitvector1', 2*(rand()-0.5), ...      % Random inclination [0, π/4]
%             'unitvector2', 2*(rand()-0.5), ...        % Random azimuth [0, 2π]
%             'unitvector3', rand(), ...        % Random azimuth [0, 2π]
%             'photons', 1000);                % Initial photon estimate
% 
%         % Fit result - now includes nested sampling results
%         estimatesPositionDefocus
%         logEvidence  % Log evidence from nested sampling
%         posteriorSamples  % Full posterior samples
%         parameterUncertainties  % Parameter uncertainties
%         credibleIntervals  % Credible intervals
% 
%         % BFP type indicator
%         model = 'hinterer'  % Default to Hinterer BackFocalPlane
% 
%         % Nested sampling parameters
%         nestedSamplingParams = struct( ...
%             'Nlive', 60, ...
%             'tolerance', 0.5, ...
%             'Nmcmc', 0, ...  % 0 for multinest, >0 for MCMC
%             'plotPosteriors', false, ...
%             'estimationMethod', 'mean')  % 'mean', 'median', or 'map'
%     end
% 
%     methods
%         function obj = FitPSF_posterior_twostage(psf, par, model, autoFit)
% 
%             if nargin < 4
%                 autoFit = true;  % Default behavior for backward compatibility
%             end
% 
%             if nargin > 2
%                 % If model is provided, set it
%                 if strcmpi(model, 'gaussian')
%                     obj.model = 'gaussian';
%                 elseif strcmpi(model, 'hinterer')
%                     obj.model = 'hinterer';
%                 elseif strcmpi(model, 'mortensen')
%                     obj.model = 'mortensen';
%                 end
%             end
% 
%             if nargin > 1
%                 obj = setInputParameters('FitPSF_posterior_twostage', obj, par);
%             end
% 
%             if nargin > 0
%                 obj.psf = psf;
%                 obj.image = psf.image;
%                 obj.nPhotonEstimate = round(sum(sum(obj.image - obj.noiseEstimate)));
%                 % Only run fitting if autoFit is true
%                 if autoFit
%                     [obj.estimatesPositionDefocus, obj.logEvidence, obj.posteriorSamples, ...
%                      obj.parameterUncertainties, obj.credibleIntervals] = fitting(obj, model);
%                 end
%             end
%         end
% 
%         %% Fit using Nested Sampling
%         function [estimatesPositionDefocus, logEvidence, posteriorSamples, parameterUncertainties, credibleIntervals] = fitting(obj, model)
%             parPsfEstimate = FitPSF_posterior_twostage.readParametersEstimate(obj.psf);
%             parPsfEstimate.dipole = Dipole(0, 0);
%             parPsfEstimate.position = Length([0 0 0], 'nm');
%             parPsfEstimate.nPhotons = obj.nPhotonEstimate;
%             parPsfEstimate.defocus = Length(0, 'nm');
%             parPsfEstimate.backgroundNoise = 0; % background noise is added later
%             parPsfEstimate.pixelSensitivityMask = obj.pixelSensitivityMask;
%             parPsfEstimate.stageDrift = obj.stageDrift;
% 
%             % Create PSF estimate based on model
%             if strcmpi(obj.model, 'hinterer')
%                 psfEstimate = PSF(parPsfEstimate);
%             elseif strcmpi(obj.model, 'mortensen')
%                 psfEstimate = PSF_mortensen(parPsfEstimate);
%             elseif strcmpi(obj.model, 'gaussian')
%                 psfEstimate = PSF_gaussian(parPsfEstimate);
%             else
%                 error('Unknown model type: %s', obj.model);
%             end
% 
%             % Prepare image
%             if strcmpi(model, 'gaussian')
%                 psfImage = obj.image ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
%             elseif strcmpi(model, 'hinterer')
%                 psfImage = obj.image;% ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
%             elseif strcmpi(model, 'mortensen')
%                 psfImage = obj.image;% ./ norm(obj.image); % dave apr 2025 - did hinterer need this?
%             end
% 
%             % Run nested sampling
%             [logEvidence, nest_samples, post_samples] = fitNestedSamplingPSF(obj, psfImage, psfEstimate, model);
% 
%             % Extract point estimates
%             estimatesPositionDefocus.ML = FitPSF_posterior_twostage.extractPointEstimates(post_samples, obj.nestedSamplingParams.estimationMethod);
% 
%             % Store full results
%             posteriorSamples = post_samples;
% 
%             % Get parameter uncertainties
%             [parameterUncertainties, credibleIntervals] = FitPSF_posterior_twostage.getParameterUncertainties(post_samples, 0.95);
%         end
% 
%         function [logZ, nest_samples, post_samples] = fitNestedSamplingPSF(obj, image, psfEstimate, model)
%             fprintf('Starting nested sampling PSF fit...\n');
% 
%             % CRITICAL FIX: Ensure ALL parameters are scalar doubles
%             Nlive = double(obj.nestedSamplingParams.Nlive);
%             if ~isscalar(Nlive)
%                 Nlive = Nlive(1);  % Take first element if it's an array
%             end
% 
%             tolerance = double(obj.nestedSamplingParams.tolerance);
%             if ~isscalar(tolerance)
%                 tolerance = tolerance(1);
%             end
% 
%             Nmcmc = double(obj.nestedSamplingParams.Nmcmc);
%             if ~isscalar(Nmcmc)
%                 Nmcmc = Nmcmc(1);
%             end
% 
%             % Additional validation to ensure they're finite scalars
%             if ~isfinite(Nlive) || Nlive < 1
%                 Nlive = 50;  % Safe default
%             end
%             if ~isfinite(tolerance) || tolerance <= 0
%                 tolerance = 0.1;  % Safe default  
%             end
%             if ~isfinite(Nmcmc) || Nmcmc < 0
%                 Nmcmc = 0;  % Safe default
%             end
% 
%             % Force to integers where needed
%             Nlive = round(Nlive);
%             Nmcmc = round(Nmcmc);
% 
%             fprintf('Validated parameters: Nlive=%d, tolerance=%.3f, Nmcmc=%d\n', Nlive, tolerance, Nmcmc);
% 
%             % % Ask user for MCMC iterations if not set
%             % if isempty(Nmcmc) || Nmcmc == 0
%             %     Nmcmc = input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
%             %     if isempty(Nmcmc)
%             %         Nmcmc = 0;
%             %     end
%             %     Nmcmc = round(double(Nmcmc));  % Ensure scalar integer
%             % end
%             Nmcmc = 0;
%             Nmcmc = round(double(Nmcmc));
% 
%             % Prepare data for nested sampling - ensure it's properly formatted
%             data = double(image);
% 
%             fprintf('Data validation:\n');
%             fprintf('Data size: %s\n', mat2str(size(data)));
%             fprintf('Data class: %s\n', class(data));
%             fprintf('Data range: [%.2e, %.2e]\n', min(data(:)), max(data(:)));
%             fprintf('Data finite: %s\n', mat2str(all(isfinite(data(:)))));
% 
%             % Define the likelihood function with additional validation
%             likelihood = @(data_in, model_func, parnames, parvals) FitPSF_posterior_twostage.validateAndCallLikelihood(obj, psfEstimate, data_in, parnames, parvals);
% 
%             % Define the PSF model function (placeholder)
%             model_func = @(data_in, parnames, parvals) ones(size(data_in));
% 
%             if strcmpi(model, 'gaussian')
%                 % Define priors for Gaussian model (4 parameters)
%                 xBounds = obj.parameterBounds.x.inNanometer;
%                 yBounds = obj.parameterBounds.y.inNanometer;
%                 defocusBounds = obj.parameterBounds.defocus.inNanometer;
%                 photonsBounds = obj.parameterBounds.photons;
% 
%                 % Ensure all bounds are scalar doubles and finite
%                 xLower = double(xBounds(1));
%                 xUpper = double(xBounds(2));
%                 yLower = double(yBounds(1));
%                 yUpper = double(yBounds(2));
%                 defocusLower = double(defocusBounds(1));
%                 defocusUpper = double(defocusBounds(2));
%                 photonsLower = double(photonsBounds(1));
%                 photonsUpper = double(photonsBounds(2));
% 
%                 % Additional safety checks
%                 if ~isfinite(xLower) || ~isfinite(xUpper) || xLower >= xUpper
%                     xLower = -100; xUpper = 100;
%                 end
%                 if ~isfinite(yLower) || ~isfinite(yUpper) || yLower >= yUpper
%                     yLower = -100; yUpper = 100;
%                 end
%                 if ~isfinite(defocusLower) || ~isfinite(defocusUpper) || defocusLower >= defocusUpper
%                     defocusLower = -10; defocusUpper = 10;
%                 end
%                 if ~isfinite(photonsLower) || ~isfinite(photonsUpper) || photonsLower >= photonsUpper
%                     photonsLower = 100; photonsUpper = 10000;
%                 end
% 
%                 prior = {'x', 'uniform', xLower, xUpper, 'fixed'; ...
%                          'y', 'uniform', yLower, yUpper, 'fixed'; ...
%                          'defocus', 'uniform', defocusLower, defocusUpper, 'fixed'; ...
%                          'photons', 'uniform', photonsLower, photonsUpper, 'fixed'};
% 
%             else
%                 % Define priors for Hinterer/Mortensen model (6 parameters: x, y, defocus, theta, phi, photons)
%                 xBounds = obj.parameterBounds.x.inNanometer;
%                 yBounds = obj.parameterBounds.y.inNanometer;
%                 defocusBounds = obj.parameterBounds.defocus.inNanometer;
%                 unitvector1Bounds = obj.parameterBounds.unitvector1;      % [0, π/2]
%                 unitvector2Bounds = obj.parameterBounds.unitvector2;      % [0, π/2]
%                 unitvector3Bounds = obj.parameterBounds.unitvector3;      % [0, π/2]
%                 photonsBounds = obj.parameterBounds.photons;
% 
%                 % Ensure all bounds are scalar doubles and finite
%                 xLower = double(xBounds(1));
%                 xUpper = double(xBounds(2));
%                 yLower = double(yBounds(1));
%                 yUpper = double(yBounds(2));
%                 defocusLower = double(defocusBounds(1));
%                 defocusUpper = double(defocusBounds(2));
%                 unitvector1Lower = double(unitvector1Bounds(1));
%                 unitvector1Upper = double(unitvector1Bounds(2));
%                 unitvector2Lower = double(unitvector2Bounds(1));
%                 unitvector2Upper = double(unitvector2Bounds(2));
%                 unitvector3Lower = double(unitvector3Bounds(1));
%                 unitvector3Upper = double(unitvector3Bounds(2));
%                 photonsLower = double(photonsBounds(1));
%                 photonsUpper = double(photonsBounds(2));
% 
%                 % Additional safety checking for all bounds
%                 bounds_check = [xLower, xUpper; yLower, yUpper; defocusLower, defocusUpper; ...
%                                unitvector1Lower, unitvector1Upper; unitvector2Lower, unitvector2Upper; photonsLower, photonsUpper];
% 
%                 for i = 1:size(bounds_check, 1)
%                     if ~isfinite(bounds_check(i,1)) || ~isfinite(bounds_check(i,2)) || bounds_check(i,1) >= bounds_check(i,2)
%                         switch i
%                             case 1  % x bounds
%                                 xLower = -100; xUpper = 100;
%                             case 2  % y bounds
%                                 yLower = -100; yUpper = 100;
%                             case 3  % defocus bounds
%                                 defocusLower = -10; defocusUpper = 10;
%                             case 4  % theta bounds
%                                 unitvector1Lower = -1; unitvector1Upper = 1;
%                             case 5  % phi bounds
%                                 unitvector2Lower = -1; unitvector2Upper = 1;
%                             case 6  % phi bounds
%                                 unitvector3Lower = 0; unitvector3Upper = 1;
%                             case 7  % photons bounds
%                                 photonsLower = 100; photonsUpper = 10000;
%                         end
%                     end
%                 end
% 
%                 prior = {'x', 'uniform', xLower, xUpper, 'fixed'; ...
%                          'y', 'uniform', yLower, yUpper, 'fixed'; ...
%                          'defocus', 'uniform', defocusLower, defocusUpper, 'fixed'; ...
%                          'unitvector1', 'uniform', unitvector1Lower, unitvector1Upper, 'fixed'; ...
%                          'unitvector2', 'uniform', unitvector2Lower, unitvector2Upper, 'fixed'; ...
%                          'unitvector3', 'uniform', unitvector3Lower, unitvector3Upper, 'fixed'; ...
%                          'photons', 'uniform', photonsLower, photonsUpper, 'fixed'};
%             end
% 
%             % Final validation of prior structure - ensure all entries are scalar
%             fprintf('Prior validation:\n');
%             for i = 1:size(prior, 1)
%                 param_name = prior{i,1};
%                 prior_type = prior{i,2};
%                 lower_bound = prior{i,3};
%                 upper_bound = prior{i,4};
% 
%                 % Force to scalar if needed
%                 if ~isscalar(lower_bound)
%                     lower_bound = lower_bound(1);
%                     prior{i,3} = lower_bound;
%                 end
%                 if ~isscalar(upper_bound)
%                     upper_bound = upper_bound(1);
%                     prior{i,4} = upper_bound;
%                 end
% 
%                 fprintf('  %s: %s [%.6f, %.6f]\n', param_name, prior_type, lower_bound, upper_bound);
% 
%                 if ~isscalar(lower_bound) || ~isscalar(upper_bound)
%                     error('Prior bounds must be scalar for parameter %s', param_name);
%                 end
%                 if ~isfinite(lower_bound) || ~isfinite(upper_bound)
%                     error('Prior bounds must be finite for parameter %s', param_name);
%                 end
%                 if lower_bound >= upper_bound
%                     error('Lower bound must be less than upper bound for parameter %s', param_name);
%                 end
%             end
% 
%             extraparams = {};  % Ensure this is empty cell array
% 
%             fprintf('\nRunning nested sampler with %d live points...\n', Nlive);
% 
%             % Run nested sampling with all scalar parameters
%             try
%                 % CRITICAL: Ensure all optional parameters are also scalar
%                 diffevfrac_val = double(2);
%                 covfrac_val = double(10);
%                 walkfrac_val = double(0);
%                 stretchfrac_val = double(0);
% 
%                 [logZ, nest_samples, post_samples] = nested_sampler(...
%                     data, ...                    % data matrix
%                     Nlive, ...                   % scalar integer  
%                     tolerance, ...               % scalar double
%                     likelihood, ...              % function handle
%                     model_func, ...              % function handle
%                     prior, ...                   % validated cell array
%                     extraparams, ...             % empty cell array
%                     'Nmcmc', Nmcmc, ...         % scalar integer
%                     'diffevfrac', diffevfrac_val, ...  % scalar double
%                     'covfrac', covfrac_val, ...        % scalar double
%                     'walkfrac', walkfrac_val, ...      % scalar double
%                     'stretchfrac', stretchfrac_val);   % scalar double
% 
%                 fprintf('Nested sampling completed successfully. Log evidence: %.2f\n', logZ);
% 
%             catch ME
%                 fprintf('\n=== NESTED SAMPLER ERROR ===\n');
%                 fprintf('Error message: %s\n', ME.message);
%                 fprintf('Error identifier: %s\n', ME.identifier);
% 
%                 for i = 1:length(ME.stack)
%                     fprintf('Stack frame %d: %s (line %d)\n', i, ME.stack(i).file, ME.stack(i).line);
%                 end
% 
%                 fprintf('============================\n');
%                 rethrow(ME);
%             end
% 
%             % Optional: Plot posterior distributions  
%             if obj.nestedSamplingParams.plotPosteriors
%                 fprintf('Plotting posterior distributions...\n');
% 
%                 try
%                     if strcmpi(model, 'gaussian')
%                         for i = 1:4
%                             posteriors(post_samples, i, {prior{i,1}});
%                         end
%                         posteriors(post_samples, [1 2], {prior{1,1}, prior{2,1}});
%                     else
%                         % For spherical coordinates (6 parameters)
%                         for i = 1:7
%                             posteriors(post_samples, i, {prior{i,1}});
%                         end
%                         posteriors(post_samples, [1 2], {prior{1,1}, prior{2,1}});
%                         % Plot angles together
%                         posteriors(post_samples, [4 5], {prior{4,1}, prior{5,1}});
%                     end
%                 catch plot_ME
%                     fprintf('Warning: Could not create posterior plots: %s\n', plot_ME.message);
%                 end
%             end
%         end
% 
%         function currentlnpdf = lnpdfFunction(obj, psfEstimate, z, lateralPositionAndDefocus) 
%             try
%                 currentPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus); 
% 
%                 % Ensure currentPSF and z have the same size
%                 if ~isequal(size(currentPSF), size(z))
%                     error('Size mismatch: currentPSF size = %s, z size = %s', ...
%                         mat2str(size(currentPSF)), mat2str(size(z)));
%                 end
% 
%                 % Compute log-likelihood element-wise, then sum
%                 log_psf = log(currentPSF);
% 
%                 % Check for invalid values
%                 if any(~isfinite(log_psf(:))) || any(~isfinite(currentPSF(:)))
%                     currentlnpdf = -inf;
%                     return;
%                 end
% 
%                 % Compute Poisson log-likelihood
%                 individual_terms = z .* log_psf - currentPSF - gammaln(z + 1);
% 
%                 % Sum all terms to get scalar result
%                 currentlnpdf = sum(individual_terms(:));
% 
%                 % Final check that result is scalar
%                 if ~isscalar(currentlnpdf)
%                     error('lnpdfFunction returned non-scalar result');
%                 end
% 
%             catch ME
%                 fprintf('Error in lnpdfFunction: %s\n', ME.message);
%                 currentlnpdf = -inf;
%             end
%         end
% 
%         function currentFitPSF = createFitPSF(obj, psfEstimate, lateralPositionAndDefocus)
% 
%             if strcmpi(obj.model, 'gaussian')
%                 % Gaussian model: x, y, defocus, photons
%                 xEstimate = lateralPositionAndDefocus(1);
%                 yEstimate = lateralPositionAndDefocus(2);
%                 photonEstimate = lateralPositionAndDefocus(4);
% 
%                 psfEstimate.position = Length([xEstimate, yEstimate, 0], 'nm');
%                 psfEstimate.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');
% 
%             else
%                 % Non-Gaussian models: x, y, defocus, theta, phi, photons
%                 xEstimate = lateralPositionAndDefocus(1);
%                 yEstimate = lateralPositionAndDefocus(2);
% 
%                 % Extract vector components
%                 unitvector1_raw = lateralPositionAndDefocus(4);
%                 unitvector2_raw = lateralPositionAndDefocus(5);
%                 unitvector3_raw = lateralPositionAndDefocus(6);
% 
%                 % Normalize to ensure unit vector
%                 norm_factor = sqrt(unitvector1_raw^2 + unitvector2_raw^2 + unitvector3_raw^2);
%                 if norm_factor < 1e-10
%                     % Handle near-zero case
%                     unitvector1_raw = 0;
%                     unitvector2_raw = 0;
%                     unitvector3_raw = 1;
%                 else
%                     unitvector1_raw = unitvector1_raw / norm_factor;
%                     unitvector2_raw = unitvector2_raw / norm_factor;
%                     unitvector3_raw = unitvector3_raw / norm_factor;
%                 end
% 
%                 % Now use normalized values to calculate inclination and azimuth
%                 inclination = acos(unitvector3_raw);
%                 azimuth = atan2(unitvector2_raw, unitvector1_raw);
%                 inclination = mod(inclination, pi/2);
%                 azimuth = mod(azimuth, 2*pi);
%                 photonEstimate = lateralPositionAndDefocus(7);
% 
%                 psfEstimate.position = Length([xEstimate, yEstimate, 0], 'nm');
%                 psfEstimate.defocus = Length(obj.parameterStartValues.defocus.inNanometer, 'nm');%Length(lateralPositionAndDefocus(3), 'nm');
%                 psfEstimate.dipole = Dipole(inclination, azimuth);
%             end
% 
%             % dave jan 2025
%             % doing more than the reduced form they were doing            
%             if strcmpi(obj.model, 'mortensen')
% 
%                 % Use Mortensen model
% 
%                 % Call Python stuff
%                 % Try two specific locations for the Python module (local vs cluster)
%                 possibleDirs = {
%                     fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
%                     '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits'
%                 };
% 
%                 % Try each location until we find the Python file
%                 pyModuleFound = false;
%                 for i = 1:length(possibleDirs)
%                     pyDir = possibleDirs{i};
% 
%                     % Check if directory exists
%                     if ~exist(pyDir, 'dir')
%                         continue;  % Skip to next directory
%                     end
% 
%                     % Check if Python file exists in this directory
%                     pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
%                     if exist(pyFilePath, 'file')
% 
%                         % Add to Python path
%                         if count(py.sys.path(), pyDir) == 0
%                             py.sys.path().insert(int32(0), pyDir);
%                         end
% 
%                         pyModuleFound = true;
%                         break;  % Found the file, stop looking
%                     end
%                 end
% 
%                 % If the module wasn't found in any location, show an error
%                 if ~pyModuleFound
%                     error(['Python module "vectorized_mortensen_flipy_phi_theta.py" not found in any of the specified locations. ', ...
%                            'Please ensure the file exists in one of these directories:\n', ...
%                            '- %s\n', ...
%                            '- %s'], possibleDirs{1}, possibleDirs{2});
%                 end
% 
%                 % Run the function defined in your python file
%                 currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
%                     xEstimate, ... % x
%                     yEstimate, ... % y
%                     inclination, ... % theta
%                     azimuth, ... % phi
%                     psfEstimate.nPixels,... % image_size_px
%                     double(psfEstimate.pixelSize.inNanometer), ... % pixel_size_nm
%                     double(psfEstimate.wavelength.inNanometer), ... % wavelength
%                     psfEstimate.refractiveIndices(2), ... % n_objective
%                     psfEstimate.refractiveIndices(1), ... % n_sample
%                     psfEstimate.objectiveNA, ... % NA
%                     photonEstimate ... % n_photons
%                 );
% 
%                 % Convert to matlab array
%                 % Get the shape of the array
%                 py_shape = currentPsf.shape;
%                 rows = double(py_shape{1});
%                 cols = double(py_shape{2});
% 
%                 % Initialize MATLAB array
%                 psf_matlab = zeros(rows, cols);
% 
%                 % Copy values individually using item() method with explicit integer conversion
%                 for i = 0:(rows-1)
%                     for j = 0:(cols-1)
%                         % Use py.int to explicitly convert indices to Python integers
%                         psf_matlab(i+1, j+1) = double(currentPsf.item(py.tuple({py.int(i), py.int(j)})));
%                     end
%                 end
% 
%                 currentPsf = psf_matlab;
% 
%                 % this bit is the equivalent of the getIntensitiesCamera()
%                 % stuff done in hinterer
%                 totalIntensity = sum(sum(currentPsf));
%                 currentPsf = currentPsf / totalIntensity * photonEstimate;
%                 currentFitPSF = currentPsf;
% 
%             else
% 
%                 % Select appropriate BackFocalPlane function based on model
%                 if strcmpi(obj.model, 'gaussian')
%                     bfp = BackFocalPlane_gaussian(psfEstimate); % Use Gaussian model
%                 elseif strcmpi(obj.model, 'hinterer')
%                     bfp = BackFocalPlane(psfEstimate); % Use Hinterer model
%                 end
% 
%                 psfEstimate.backFocalPlane = bfp;
% 
%                 % dave
%                 psfEstimate.fieldBFP.x = bfp.electricField.x;
%                 psfEstimate.fieldBFP.y = bfp.electricField.y;
% 
%                 currentPsf = zeros(psfEstimate.nPixels,psfEstimate.nPixels); 
%                 for k=1:size(psfEstimate.stageDrift.motion,1)
%                     % Apply aberrations
%                     aberrationCoeffs = getAberrations(psfEstimate,k);
%                     fieldBFP = applyAberrations(psfEstimate, aberrationCoeffs);
%                     % Get image from BFP field
%                     currentPsf = currentPsf + getIntensitiesCamera(psfEstimate, fieldBFP);
%                 end
% 
%                 totalIntensity = sum(currentPsf,'all');
%                 currentPsf = currentPsf ./ totalIntensity * photonEstimate + obj.noiseEstimate;
%                 currentFitPSF = currentPsf;% !!! maybe swap this back - 8th May 2025
% 
%             end
% 
%             currentFitPSF = adjustExcitation(psfEstimate, currentFitPSF);
%             currentFitPSF = applyShotNoise(psfEstimate, currentFitPSF);
%             currentFitPSF = addBackgroundNoise(psfEstimate, currentFitPSF);
%         end
% 
%         % Method to set nested sampling parameters
%         function obj = setNestedSamplingParams(obj, varargin)
%             % Set nested sampling parameters
%             % Usage: obj = obj.setNestedSamplingParams('Nlive', 1000, 'tolerance', 0.01, ...)
% 
%             for i = 1:2:length(varargin)
%                 param = varargin{i};
%                 value = varargin{i+1};
% 
%                 if isfield(obj.nestedSamplingParams, param)
%                     obj.nestedSamplingParams.(param) = value;
%                 else
%                     warning('Unknown nested sampling parameter: %s', param);
%                 end
%             end
%         end
% 
%         function obj = runFittingWithParams(obj, varargin)
%             % Set parameters first
%             obj = obj.setNestedSamplingParams(varargin{:});
% 
%             % Then run fitting
%             [obj.estimatesPositionDefocus, obj.logEvidence, obj.posteriorSamples, ...
%              obj.parameterUncertainties, obj.credibleIntervals] = fitting(obj, obj.model);
%         end
% 
%         % Method to get parameter summary
%         function summary = getParameterSummary(obj)
%             % Get a summary of parameter estimates and uncertainties
% 
%             if isempty(obj.posteriorSamples)
%                 error('No posterior samples available. Run fitting first.');
%             end
% 
%             samples = obj.posteriorSamples.samples;
% 
%             if strcmpi(obj.model, 'gaussian')
%                 param_names = {'x', 'y', 'defocus', 'photons'};
%             else
%                 param_names = {'x', 'y', 'defocus', 'unitvector1', 'unitvector2', 'unitvector3', 'photons'};
%             end
% 
%             summary = struct();
%             summary.parameter_names = param_names;
%             summary.estimates = obj.estimatesPositionDefocus.ML(1:end-1);  % Exclude cost
%             summary.uncertainties = obj.parameterUncertainties;
%             summary.credible_intervals = obj.credibleIntervals;
%             summary.log_evidence = obj.logEvidence;
%             summary.estimation_method = obj.nestedSamplingParams.estimationMethod;
% 
%             % Display summary
%             fprintf('\n=== Parameter Estimation Summary ===\n');
%             fprintf('Model: %s\n', obj.model);
%             fprintf('Estimation method: %s\n', summary.estimation_method);
%             fprintf('Log evidence: %.2f\n', summary.log_evidence);
%             fprintf('\nParameter estimates (±1σ, [95%% CI]):\n');
% 
%             for i = 1:length(param_names)
%                 fprintf('%s: %.2f ± %.2f [%.2f, %.2f]\n', ...
%                     param_names{i}, ...
%                     summary.estimates(i), ...
%                     summary.uncertainties(i), ...
%                     summary.credible_intervals(i,1), ...
%                     summary.credible_intervals(i,2));
%             end
%             fprintf('\n');
%         end
%     end
% 
%     methods (Static)
%         % function par = readParametersEstimate(psf)
%         %     % dave apr 2025 - add this so we can inherit the right params when
%         %     % fitting
%         %     par = struct();
%         %     par.nPixels = psf.nPixels;
%         %     par.pixelSize = psf.pixelSize;
%         %     par.wavelength = psf.wavelength;
%         %     par.objectiveNA = psf.objectiveNA;
%         %     par.pixelSensitivityMask = psf.pixelSensitivityMask;
%         %     par.nDiscretizationBFP = psf.nDiscretizationBFP;
%         % 
%         %     % Fluorophore
%         %     par.shotNoise = 0; 
%         %     par.reducedExcitation = 0;
%         % 
%         %     % Microscope setup
%         %     par.wavelength = psf.wavelength;
%         %     if isprop(psf, 'astigmatism') || isfield(psf, 'astigmatism')
%         %         par.astigmatism = psf.astigmatism;
%         %     end
%         %     par.objectiveNA = psf.objectiveNA;
%         %     par.objectiveFocalLength = psf.objectiveFocalLength;
%         %     par.refractiveIndices = psf.refractiveIndices;
%         %     if isprop(psf, 'heightIntermediateLayer') || isfield(psf, 'heightIntermediateLayer')
%         %         par.heightIntermediateLayer = psf.heightIntermediateLayer;
%         %     end
%         % 
%         %     % Back focal plane
%         %     if isprop(psf, 'phaseMask') || isfield(psf, 'phaseMask')
%         %         par.phaseMask = psf.phaseMask;
%         %     end
%         %     par.nDiscretizationBFP = psf.nDiscretizationBFP;
%         % 
%         %     % Camera
%         %     par.pixelSize = psf.pixelSize;
%         % 
%         %     % dave apr 2025
%         %     par.nPixels = psf.nPixels;
%         % end
% 
%         function logL = validateAndCallLikelihood(obj, psfEstimate, data, parnames, parvals)
%             try
%                 % Simple validation
%                 if length(parnames) ~= length(parvals)
%                     logL = -inf;
%                     return;
%                 end
% 
%                 % Call the nested likelihood function - note the argument order
%                 logL = FitPSF_posterior_twostage.lnpdfFunction_nested(obj, obj, psfEstimate, data, parnames, parvals);
% 
%                 % Final validation
%                 if ~isscalar(logL) || ~isreal(logL)
%                     fprintf('Warning: Likelihood function returned invalid result\n');
%                     logL = -inf;
%                 end
% 
%             catch ME
%                 fprintf('Error in validateAndCallLikelihood: %s\n', ME.message);
%                 if length(ME.stack) > 0
%                     fprintf('  at line %d in %s\n', ME.stack(1).line, ME.stack(1).file);
%                 end
%                 logL = -inf;
%             end
%         end
% 
%         function logL = lnpdfFunction_nested(obj, obj_ref, psfEstimate, data, parnames, parvals)
%             % Convert parameter names and values to the format expected by your original function
% 
%             % Debug: Print what we're receiving
%             persistent debug_count;
%             if isempty(debug_count)
%                 debug_count = 0;
%             end
%             debug_count = debug_count + 1;
% 
%             try
%                 % Handle the case where nested_sampler might pass arguments differently
%                 % Sometimes nested_sampler doesn't follow the expected format exactly
% 
%                 % If data is not a matrix, it might be that the arguments are in a different order
%                 if ~ismatrix(data) || ~isnumeric(data)
%                     % Try to find the actual data in the arguments
%                     % Sometimes the first argument might be something else
%                     fprintf('Warning: data is not a numeric matrix, trying to find actual data...\n');
% 
%                     % Check if obj_ref contains the data we need
%                     if isfield(obj_ref, 'image') && isnumeric(obj_ref.image)
%                         actual_data = obj_ref.image;
%                         fprintf('Found data in obj_ref.image: %s\n', mat2str(size(actual_data)));
%                     elseif isnumeric(obj_ref) && ismatrix(obj_ref)
%                         actual_data = obj_ref;
%                         fprintf('Using obj_ref as data: %s\n', mat2str(size(actual_data)));
%                     else
%                         % Last resort: use the original data from obj
%                         actual_data = obj.image;
%                         fprintf('Using obj.image as fallback: %s\n', mat2str(size(actual_data)));
%                     end
%                 else
%                     actual_data = data;
%                 end
% 
%                 % Ensure we have valid data
%                 if ~ismatrix(actual_data) || ~isnumeric(actual_data)
%                     logL = -inf;
%                     return;
%                 end
% 
%                 % Build parameter vector in the expected order
%                 if length(parnames) == 4  % Gaussian model
%                     params = zeros(1, 4);
%                     for i = 1:length(parnames)
%                         val = parvals{i};
%                         if ~isscalar(val) || ~isreal(val) || ~isfinite(val)
%                             logL = -inf;
%                             return;
%                         end
% 
%                         switch parnames{i}
%                             case 'x'
%                                 params(1) = double(val);
%                             case 'y'
%                                 params(2) = double(val);
%                             case 'defocus'
%                                 params(3) = double(val);
%                             case 'photons'
%                                 params(4) = double(val);
%                         end
%                     end
%                 else  % Full model with spherical angles (6 parameters)
%                     params = zeros(1, 7);
%                     for i = 1:length(parnames)
%                         val = parvals{i};
%                         if ~isscalar(val) || ~isreal(val) || ~isfinite(val)
%                             logL = -inf;
%                             return;
%                         end
% 
%                         switch parnames{i}
%                             case 'x'
%                                 params(1) = double(val);
%                             case 'y'
%                                 params(2) = double(val);
%                             case 'defocus'
%                                 params(3) = double(val);
%                             case 'unitvector1'
%                                 params(4) = double(val);
%                             case 'unitvector2'
%                                 params(5) = double(val);
%                             case 'unitvector3'
%                                 params(6) = double(val);
%                             case 'photons'
%                                 params(7) = double(val);
%                         end
%                     end
%                 end
% 
%                 % Additional parameter validation
%                 if length(params) == 4 && params(4) <= 0  % Gaussian model photons
%                     logL = -inf;
%                     return;
%                 elseif length(params) == 6 && params(6) <= 0  % Full model photons
%                     logL = -inf;
%                     return;
%                 end
% 
%                 % Call your original likelihood function with the correct object reference
%                 % Use obj (the FitPSF_posterior_twostage object) instead of obj_ref
%                 logL_result = obj.lnpdfFunction(psfEstimate, actual_data, params);
% 
%                 % CRITICAL: Ensure the result is a scalar
%                 if ~isscalar(logL_result)
%                     error('Likelihood function returned non-scalar result: size = %s', mat2str(size(logL_result)));
%                 end
% 
%                 % Ensure log-likelihood is finite and real
%                 if ~isreal(logL_result) || ~isfinite(logL_result)
%                     logL = -inf;
%                 else
%                     logL = double(logL_result);  % Ensure it's a double scalar
%                 end
% 
%             catch ME
%                 % If any error occurs in likelihood computation, return very low likelihood
%                 fprintf('Error in likelihood computation: %s\n', ME.message);
%                 if length(ME.stack) > 0
%                     fprintf('Error in file %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
%                 end
%                 logL = -inf;
%             end
%         end
% 
%         function model_output = psfModel_nested(obj, obj_ref, psfEstimate, data, parnames, parvals)
%             % This function isn't strictly necessary for your PSF fitting since you're
%             % computing likelihood directly from the image, but nested_sampler requires it
%             model_output = ones(size(data));  % Placeholder
%         end
% 
%         function estimatesPositionDefocusML = extractPointEstimates(post_samples, method)
%             % Extract point estimates from posterior samples
%             % 
%             % Inputs:
%             %   post_samples - posterior samples from nested sampling
%             %   method - 'mean', 'median', 'map' (maximum a posteriori)
%             %
%             % Output:
%             %   estimatesPositionDefocusML - parameter estimates in same format as original
% 
%             if nargin < 2
%                 method = 'mean';  % Default to posterior mean
%             end
% 
%             % Handle different formats of post_samples
%             if isstruct(post_samples) && isfield(post_samples, 'samples') && isfield(post_samples, 'logL')
%                 samples = post_samples.samples;
%                 logL_values = post_samples.logL;
%             elseif iscell(post_samples) && length(post_samples) >= 2
%                 % Alternative format from nested_sampler
%                 samples = post_samples{1};
%                 logL_values = post_samples{2};
%             else
%                 % Try to adapt to whatever format we have
%                 if isnumeric(post_samples)
%                     % Assume post_samples is just the samples matrix and we don't have logL values
%                     samples = post_samples;
%                     logL_values = zeros(size(samples, 1), 1);  % Dummy values
%                     fprintf('Warning: No log-likelihood values available. Using dummy values.\n');
%                 else
%                     error('Unrecognized format for post_samples');
%                 end
%             end
% 
%             switch lower(method)
%                 case 'mean'
%                     % Posterior mean (most common choice)
%                     estimatesPositionDefocusML = mean(samples, 1);
% 
%                 case 'median'
%                     % Posterior median (more robust to outliers)
%                     estimatesPositionDefocusML = median(samples, 1);
% 
%                 case 'map'
%                     % Maximum a posteriori (MAP) estimate - sample with highest likelihood
%                     [~, map_idx] = max(logL_values);
%                     estimatesPositionDefocusML = samples(map_idx, :);
% 
%                 otherwise
%                     error('Unknown method. Use ''mean'', ''median'', or ''map''');
%             end
% 
%             % Add the log-likelihood value at the end (equivalent to cost in original)
%             % For MAP estimate, use the actual max likelihood
%             % For mean/median, evaluate likelihood at that point
%             if strcmpi(method, 'map') && ~isempty(logL_values)
%                 estimatesPositionDefocusML(end+1) = -logL_values(map_idx);  % Convert to cost (negative log-likelihood)
%             else
%                 % For mean/median, we use the mean of the log-likelihood values
%                 if ~isempty(logL_values)
%                     estimatesPositionDefocusML(end+1) = -mean(logL_values);
%                 else
%                     estimatesPositionDefocusML(end+1) = 0;  % Default value if no log-likelihood available
%                 end
%             end
%         end
% 
%         function [param_uncertainties, credible_intervals] = getParameterUncertainties(post_samples, confidence_level)
%             % Extract parameter uncertainties from posterior samples
%             %
%             % Inputs:
%             %   post_samples - posterior samples from nested sampling  
%             %   confidence_level - confidence level for credible intervals (default: 0.95)
%             %
%             % Outputs:
%             %   param_uncertainties - standard deviations of parameters
%             %   credible_intervals - credible intervals [lower, upper] for each parameter
% 
%             if nargin < 2
%                 confidence_level = 0.95;
%             end
% 
%             % Handle different formats of post_samples
%             if isstruct(post_samples) && isfield(post_samples, 'samples')
%                 samples = post_samples.samples;
%             elseif iscell(post_samples) && length(post_samples) >= 1
%                 % Alternative format from nested_sampler
%                 samples = post_samples{1};
%             else
%                 % Try to adapt to whatever format we have
%                 if isnumeric(post_samples)
%                     % Assume post_samples is just the samples matrix
%                     samples = post_samples;
%                 else
%                     error('Unrecognized format for post_samples');
%                 end
%             end
% 
%             % Parameter uncertainties (standard deviations)
%             param_uncertainties = std(samples, 1);
% 
%             % Credible intervals
%             alpha = 1 - confidence_level;
%             lower_percentile = 100 * alpha / 2;
%             upper_percentile = 100 * (1 - alpha / 2);
% 
%             credible_intervals = [prctile(samples, lower_percentile, 1); ...
%                                  prctile(samples, upper_percentile, 1)]';
%         end
%     end
% end