function [covarMatrix, fisherMatrix] = calculateSphericalCoordinatesCovarianceMatrix(paramEst, x_est, y_est, defocus_est, theta_est, phi_est, photon_est, noise_est, model)
    % Calculates the theoretical covariance matrix for parameter estimates in spherical coordinates
    % based on the Fisher Information Matrix
    %
    % INPUTS:
    %   paramEst - PSF object with image and parameters
    %   x_est, y_est - Position estimates
    %   defocus_est - Defocus estimate
    %   theta_est, phi_est - Spherical angle estimates
    %   photon_est, noise_est - Photon and noise estimates
    %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
    %
    % OUTPUTS:
    %   covarMatrix - Covariance matrix (inverse of Fisher matrix)
    %   fisherMatrix - Fisher information matrix

    % Calculate probability distribution (PSF model)
    psf = calculateProbabilitySpherical(paramEst, x_est, y_est, defocus_est, theta_est, phi_est, photon_est, noise_est, model);

    % Calculate derivatives of PSF with respect to parameters
    derivatives = calculateDerivativesSpherical(paramEst, x_est, y_est, defocus_est, theta_est, phi_est, photon_est, noise_est, model);

    % Construct Fisher Information Matrix
    fisherMatrix = zeros(length(derivatives), length(derivatives));

    % Fill Fisher Matrix
    for i = 1:length(derivatives)
        for j = 1:i  % Use symmetry to save calculations
            % Formula for Fisher Matrix: I_ij = N * sum( (∂p/∂θ_i * ∂p/∂θ_j) / (p + b/N) )
            % Where p is the probability, b is background, N is photon count
            denom = psf + noise_est/photon_est;
            fisherMatrix(i,j) = sum(sum(derivatives{i} .* derivatives{j} ./ denom));

            % Matrix is symmetric
            if i ~= j
                fisherMatrix(j,i) = fisherMatrix(i,j);
            end
        end
    end

    % Calculate covariance matrix (Cramér-Rao lower bound)
    covarMatrix = inv(fisherMatrix);
end

function currentFitPSF = calculateProbabilitySpherical(paramEst, x, y, defocus, theta, phi, photons, noise_est, model)
    % % This function is like calculateProbability, but takes theta and phi directly
    % % Convert spherical to Cartesian representation
    % theta = mod(theta, pi);         % Ensure theta is in [0, pi]
    % phi = mod(phi, 2*pi);           % Ensure phi is in [0, 2*pi]
    % 
    % newangle1 = sin(theta) * cos(phi);
    % newangle2 = sin(theta) * sin(phi);
    % newangle3 = cos(theta);
    % 
    % % Now call the existing function with these Cartesian coordinates
    % currentFitPSF = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);

    if strcmpi(model, 'gaussian')

        paramEst.position = Length([x, y, 0], 'nm');
        paramEst.defocus = Length(defocus, 'nm');

    else

        paramEst.position = Length([x, y, 0], 'nm');
        paramEst.defocus = Length(defocus, 'nm');
        paramEst.dipole = Dipole(theta, phi);

    end

    if strcmpi(model, 'mortensen')

        % Use Mortensen model

        % Call Python stuff
        % pyDir = pwd; % or specify the exact path where mortensen_simulator.py is located
%        pyDir = '/home/tfq96423/Documents/cryoCLEM/dipole-issue/fixed-dipole-issue/hinterer/simulation-and-fit/@FitPSF_ML_reparam2'; % or specify the exact path where mortensen_simulator.py is located
        % pyDir = fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits');
        % if count(py.sys.path(), pyDir) == 0
        %     py.sys.path().insert(int32(0), pyDir);
        % end




        % Try two specific locations for the Python module (local vs
        % cluster)
        possibleDirs = {
            fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
            '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits' 
        };

        % Try each location until we find the Python file
        pyModuleFound = false;
        for i = 1:length(possibleDirs)
            pyDir = possibleDirs{i};
            % fprintf('Checking for Python module in: %s\n', pyDir);

            % Check if directory exists
            if ~exist(pyDir, 'dir')
                % fprintf('Directory does not exist: %s\n', pyDir);
                continue;  % Skip to next directory
            end

            % Check if Python file exists in this directory
            pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
            if exist(pyFilePath, 'file')
                % fprintf('Found Python module at: %s\n', pyFilePath);

                % Add to Python path
                if count(py.sys.path(), pyDir) == 0
                    py.sys.path().insert(int32(0), pyDir);
                    % fprintf('Added to Python path: %s\n', pyDir);
                end

                pyModuleFound = true;
                break;  % Found the file, stop looking
            else
                % fprintf('Python file not found in directory: %s\n', pyDir);
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
        % currentPsf = py.mortensen_simulator.run_simulator( ...
        currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
            x, ... % x
            y, ... % y
            theta, ... % theta
            phi, ... % phi
            paramEst.nPixels,... % image_size_px
            double(paramEst.pixelSize.inNanometer), ... % pixel_size_nm
            double(paramEst.wavelength.inNanometer), ... % wavelength
            paramEst.refractiveIndices(2), ... % n_objective
            paramEst.refractiveIndices(1), ... % n_sample
            paramEst.objectiveNA, ... % NA
            photons ... % n_photons
        );

        % Convert to matlab array
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
        currentPsf = currentPsf / totalIntensity * photons;
        currentFitPSF = currentPsf;

    else

        % Select appropriate BackFocalPlane function based on model
        if strcmpi(model, 'gaussian')
            bfp = BackFocalPlane_gaussian(paramEst); % Use Gaussian model
        elseif strcmpi(model, 'hinterer')
            bfp = BackFocalPlane(paramEst); % Use Hinterer model
        end

            paramEst.backFocalPlane = bfp;

            % dave
            paramEst.fieldBFP.x = bfp.electricField.x;
            paramEst.fieldBFP.y = bfp.electricField.y;

            currentPsf = zeros(paramEst.nPixels,paramEst.nPixels); 
            for k=1:size(paramEst.stageDrift.motion,1)
                % Apply aberrations
                aberrationCoeffs = getAberrations(paramEst,k);
                fieldBFP = applyAberrations(paramEst, aberrationCoeffs);
                % Get image from BFP field
                currentPsf = currentPsf + getIntensitiesCamera(paramEst, fieldBFP);
            end

            totalIntensity = sum(currentPsf,'all');
            % disp(totalIntensity)
            currentPsf = currentPsf ./ totalIntensity * photons + noise_est;
            % disp(photons)
            % disp(currentPsf)
            % currentFitPSF = currentPsf ./ norm(currentPsf);
            currentFitPSF = currentPsf;
            % disp(currentFitPSF)

    end

    currentFitPSF = adjustExcitation(paramEst, currentFitPSF);
    currentFitPSF = applyShotNoise(paramEst, currentFitPSF);
    currentFitPSF = addBackgroundNoise(paramEst, currentFitPSF);




end

function derivatives = calculateDerivativesSpherical(paramEst, x, y, defocus, theta, phi, photons, noise_est, model)
    % Calculate derivatives of the PSF with respect to each parameter including theta and phi
    % Using finite differences for approximation

    % Small step sizes for numerical derivatives
    delta_pos = 1;     % nm
    delta_theta = 0.01; % radians
    delta_phi = 0.01;   % radians

    % Calculate PSF at current parameters
    psf0 = calculateProbabilitySpherical(paramEst, x, y, defocus, theta, phi, photons, noise_est, model);

    % Derivatives with respect to x and y position
    psf_dx_plus = calculateProbabilitySpherical(paramEst, x + delta_pos, y, defocus, theta, phi, photons, noise_est, model);
    psf_dx_minus = calculateProbabilitySpherical(paramEst, x - delta_pos, y, defocus, theta, phi, photons, noise_est, model);
    dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);

    psf_dy_plus = calculateProbabilitySpherical(paramEst, x, y + delta_pos, defocus, theta, phi, photons, noise_est, model);
    psf_dy_minus = calculateProbabilitySpherical(paramEst, x, y - delta_pos, defocus, theta, phi, photons, noise_est, model);
    dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);

    % Derivative with respect to theta (with proper wrapping)
    % If close to boundaries (0 or pi), use one-sided differences
    if theta < delta_theta
        % Near lower boundary, use forward difference
        psf_dtheta_plus = calculateProbabilitySpherical(paramEst, x, y, defocus, theta + delta_theta, phi, photons, noise_est, model);
        dtheta = (psf_dtheta_plus - psf0) / delta_theta;
    elseif theta > (pi/2) - delta_theta
        % Near upper boundary, use backward difference
        psf_dtheta_minus = calculateProbabilitySpherical(paramEst, x, y, defocus, theta - delta_theta, phi, photons, noise_est, model);
        dtheta = (psf0 - psf_dtheta_minus) / delta_theta;
    else
        % Use central difference
        psf_dtheta_plus = calculateProbabilitySpherical(paramEst, x, y, defocus, theta + delta_theta, phi, photons, noise_est, model);
        psf_dtheta_minus = calculateProbabilitySpherical(paramEst, x, y, defocus, theta - delta_theta, phi, photons, noise_est, model);
        dtheta = (psf_dtheta_plus - psf_dtheta_minus) / (2 * delta_theta);
    end

    % Derivative with respect to phi (with proper wrapping)
    % phi is periodic with period 2*pi
    psf_dphi_plus = calculateProbabilitySpherical(paramEst, x, y, defocus, theta, mod(phi + delta_phi, 2*pi), photons, noise_est, model);
    psf_dphi_minus = calculateProbabilitySpherical(paramEst, x, y, defocus, theta, mod(phi - delta_phi, 2*pi), photons, noise_est, model);
    dphi = (psf_dphi_plus - psf_dphi_minus) / (2 * delta_phi);

    % For models that only use position (like Gaussian)
    if strcmpi(model, 'gaussian')
        derivatives = {dx, dy};
    else
        % Include orientation derivatives for other models
        derivatives = {dx, dy, dtheta, dphi};
    end
end