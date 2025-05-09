function [covarMatrix, fisherMatrix] = calculateCovarianceMatrix(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model)
    % Calculates the theoretical covariance matrix for parameter estimates
    % based on the Fisher Information Matrix
    %
    % INPUTS:
    %   psfInit - PSF object with image and parameters
    %   fitResult - Results from the PSF fitting
    %   paramEst - Parameter estimates and bounds
    %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
    %
    % OUTPUTS:
    %   covarMatrix - Covariance matrix (inverse of Fisher matrix)
    %   fisherMatrix - Fisher information matrix

    % Calculate probability distribution (PSF model)
    psf = calculateProbability(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);

    % Calculate derivatives of PSF with respect to parameters
    derivatives = calculateDerivatives(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);

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

    % % Scale by photon count
    % fisherMatrix = photon_est * fisherMatrix;

    % Calculate covariance matrix (Cramér-Rao lower bound)
    covarMatrix = inv(fisherMatrix);
end

function currentFitPSF = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)

    % Validate parameters to ensure they're in valid ranges
    newangle1 = max(-1, min(1, newangle1));
    newangle2 = max(-1, min(1, newangle2));
    newangle3 = max(0, min(1, newangle3));

    if strcmpi(model, 'gaussian')

        paramEst.position = Length([x, y, 0], 'nm');
        paramEst.defocus = Length(defocus, 'nm');

    else

        paramEst.position = Length([x, y, 0], 'nm');
        paramEst.defocus = Length(defocus, 'nm');
        theta = acos(newangle3);
        phi = atan2(newangle2, newangle1);
        theta = mod(theta, pi/2);
        phi = mod(phi, 2*pi);
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

function derivatives = calculateDerivatives(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
    % Calculate derivatives of the PSF with respect to each parameter
    % Using finite differences for approximation

    % Small step sizes for numerical derivatives
    delta_pos = 1; % nm
    delta_newangle = 0.01;
    % delta_defocus = 1; % nm
    % delta_photon = max(5, 0.05 * photons); % 5% of photon count or at least 5

    % Calculate PSF at current parameters
    psf0 = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);

    % Derivatives with respect to x and y position
    psf_dx_plus = calculateProbability(paramEst, x + delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    psf_dx_minus = calculateProbability(paramEst, x - delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);

    psf_dy_plus = calculateProbability(paramEst, x, y + delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    psf_dy_minus = calculateProbability(paramEst, x, y - delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);

    % % Derivative with respect to defocus
    % psf_dz_plus = calculateProbability(paramEst, x, y, defocus + delta_defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    % psf_dz_minus = calculateProbability(paramEst, x, y, defocus - delta_defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    % dz = (psf_dz_plus - psf_dz_minus) / (2 * delta_defocus);

    % % Derivative with respect to photon count
    % psf_dphoton_plus = calculateProbability(paramEst, x, y, defocus, theta, phi, photons + delta_photon, noise_est, model);
    % psf_dphoton_minus = calculateProbability(paramEst, x, y, defocus, theta, phi, photons - delta_photon, noise_est, model);
    % psf_dphoton_plus = psf_dphoton_plus / sum(sum(psf_dphoton_plus)) * (photons + delta_photon);
    % psf_dphoton_minus = psf_dphoton_minus / sum(sum(psf_dphoton_minus)) * (photons - delta_photon);
    % dphoton = (psf_dphoton_plus - psf_dphoton_minus) / (2 * delta_photon);

    if strcmpi(model, 'gaussian')
        % For Gaussian model, return only position and defocus derivatives
        derivatives = {dx, dy};%{dx, dy, dz, dphoton};
    else

        % Also need to do orientations if using Hinterer or Mortensen

        % For newangle1 (range [-1,1])
        if newangle1 < -1 + delta_newangle
            psf_dnewangle1_plus = calculateProbability(paramEst, x, y, defocus, newangle1 + delta_newangle, newangle2, newangle3, photons, noise_est, model);
            dnewangle1 = (psf_dnewangle1_plus - psf0) / delta_newangle;
        elseif newangle1 > 1 - delta_newangle
            psf_dnewangle1_minus = calculateProbability(paramEst, x, y, defocus, newangle1 - delta_newangle, newangle2, newangle3, photons, noise_est, model);
            dnewangle1 = (psf0 - psf_dnewangle1_minus) / delta_newangle;
        else
            psf_dnewangle1_plus = calculateProbability(paramEst, x, y, defocus, newangle1 + delta_newangle, newangle2, newangle3, photons, noise_est, model);
            psf_dnewangle1_minus = calculateProbability(paramEst, x, y, defocus, newangle1 - delta_newangle, newangle2, newangle3, photons, noise_est, model);
            dnewangle1 = (psf_dnewangle1_plus - psf_dnewangle1_minus) / (2 * delta_newangle);
        end

        % For newangle2 (range [-1,1])
        if newangle2 < -1 + delta_newangle
            psf_dnewangle2_plus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2 + delta_newangle, newangle3, photons, noise_est, model);
            dnewangle2 = (psf_dnewangle2_plus - psf0) / delta_newangle;
        elseif newangle2 > 1 - delta_newangle
            psf_dnewangle2_minus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2 - delta_newangle, newangle3, photons, noise_est, model);
            dnewangle2 = (psf0 - psf_dnewangle2_minus) / delta_newangle;
        else
            psf_dnewangle2_plus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2 + delta_newangle, newangle3, photons, noise_est, model);
            psf_dnewangle2_minus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2 - delta_newangle, newangle3, photons, noise_est, model);
            dnewangle2 = (psf_dnewangle2_plus - psf_dnewangle2_minus) / (2 * delta_newangle);
        end

        % For newangle3 (range [0,1])
        if newangle3 < delta_newangle
            % Too close to lower boundary, use forward difference only
            psf_dnewangle3_plus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_newangle, photons, noise_est, model);
            dnewangle3 = (psf_dnewangle3_plus - psf0) / delta_newangle;
        elseif newangle3 > 1 - delta_newangle
            % Too close to upper boundary, use backward difference only
            psf_dnewangle3_minus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_newangle, photons, noise_est, model);
            dnewangle3 = (psf0 - psf_dnewangle3_minus) / delta_newangle;
        else
            % Safe to use central difference
            psf_dnewangle3_plus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_newangle, photons, noise_est, model);
            psf_dnewangle3_minus = calculateProbability(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_newangle, photons, noise_est, model);
            dnewangle3 = (psf_dnewangle3_plus - psf_dnewangle3_minus) / (2 * delta_newangle);
        end

        derivatives = {dx, dy, dnewangle1, dnewangle2, dnewangle3};
    end

end