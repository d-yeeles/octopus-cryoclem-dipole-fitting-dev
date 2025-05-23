% function [cartesianCovarMatrix, sphericalCovarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model)
%     % Calculates the theoretical covariance matrix for parameter estimates in Cartesian coordinates
%     % based on the Fisher Information Matrix, and optionally converts to spherical coordinates
%     %
%     % INPUTS:
%     %   paramEst - PSF object with image and parameters
%     %   x_est, y_est - Position estimates
%     %   defocus_est - Defocus estimate
%     %   newangle1_est, newangle2_est, newangle3_est - Cartesian dipole orientation estimates
%     %       newangle1 ∈ [-1,1], newangle2 ∈ [-1,1], newangle3 ∈ [0,1]
%     %   photon_est, noise_est - Photon and noise estimates
%     %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
%     %
%     % OUTPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (inverse of Fisher matrix)
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (theta, phi)
%     %   fisherMatrix - Fisher information matrix
% 
%     % Calculate probability distribution (PSF model)
%     psf = calculateProbabilityCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Calculate derivatives of PSF with respect to parameters
%     derivatives = calculateDerivativesCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Construct Fisher Information Matrix
%     fisherMatrix = zeros(length(derivatives), length(derivatives));
% 
%     % Fill Fisher Matrix
%     for i = 1:length(derivatives)
%         for j = 1:i  % Use symmetry to save calculations
%             % Formula for Fisher Matrix: I_ij = N * sum( (∂p/∂θ_i * ∂p/∂θ_j) / (p + b/N) )
%             % Where p is the probability, b is background, N is photon count
%             denom = psf + noise_est/photon_est;
%             fisherMatrix(i,j) = sum(sum(derivatives{i} .* derivatives{j} ./ denom));
% 
%             % Matrix is symmetric
%             if i ~= j
%                 fisherMatrix(j,i) = fisherMatrix(i,j);
%             end
%         end
%     end
% 
%     % Calculate Cartesian covariance matrix (Cramér-Rao lower bound)
%     cartesianCovarMatrix = inv(fisherMatrix);
% 
%     % Convert to spherical coordinates (if not using Gaussian model)
%     if ~strcmpi(model, 'gaussian') && nargout > 1
%         sphericalCovarMatrix = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1_est, newangle2_est, newangle3_est);
%     else
%         % For Gaussian model, just return the same covariance
%         sphericalCovarMatrix = cartesianCovarMatrix;
%     end
% end
% 
% function [sphericalCovarMatrix] = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1, newangle2, newangle3)
%     % Converts a covariance matrix from Cartesian coordinates (newangle1, newangle2, newangle3)
%     % to spherical coordinates (theta, phi)
%     %
%     % INPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (5x5 or larger)
%     %                          Order is assumed to be [x, y, newangle1, newangle2, newangle3, ...]
%     %   newangle1, newangle2, newangle3 - Current values of Cartesian coordinates
%     %
%     % OUTPUTS:
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (4x4 or larger)
%     %                          Order is [x, y, theta, phi, ...]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize to ensure we have a unit vector
%     newangle1_norm = newangle1 / normFactor;
%     newangle2_norm = newangle2 / normFactor;
%     newangle3_norm = newangle3 / normFactor;
% 
%     % Calculate current theta and phi from the Cartesian coordinates
%     theta = acos(newangle3_norm);
%     phi = atan2(newangle2_norm, newangle1_norm);
% 
%     % Ensure phi is in [0, 2π]
%     if phi < 0
%         phi = phi + 2*pi;
%     end
% 
%     % Calculate the Jacobian matrix for the transformation
%     % J_ij = ∂(spherical_i)/∂(cartesian_j)
%     J = calculateJacobian(newangle1, newangle2, newangle3);
% 
%     % Extract the orientation part of the covariance matrix (3x3 submatrix)
%     % This assumes the order in cartesianCovarMatrix is [x, y, newangle1, newangle2, newangle3, ...]
%     orientationCovar = cartesianCovarMatrix(3:5, 3:5);
% 
%     % Transform the orientation covariance using the Jacobian
%     % Cov_spherical = J * Cov_cartesian * J'
%     newOrientationCovar = J * orientationCovar * J';
% 
%     % Create the new spherical covariance matrix, preserving other elements
%     sphericalCovarMatrix = cartesianCovarMatrix;
% 
%     % Replace the orientation part (3:5, 3:5) with the transformed 2x2 matrix (theta, phi)
%     % We'll need to adjust the matrix size since we're going from 3 parameters to 2
% 
%     % First, handle the case where we have a 5x5 or larger matrix
%     n = size(cartesianCovarMatrix, 1);
%     if n >= 5
%         % Create a new matrix of appropriate size (n-1 x n-1)
%         temp = zeros(n-1, n-1);
% 
%         % Copy x, y covariances directly
%         temp(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
% 
%         % Insert the transformed theta, phi covariances
%         temp(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
% 
%         % Copy cross-covariances between x,y and orientation parameters
%         % For x,y vs theta,phi
%         temp(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
%         temp(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
% 
%         % Copy any remaining parameters (if matrix is larger than 5x5)
%         if n > 5
%             temp(1:2, 5:n-1) = cartesianCovarMatrix(1:2, 6:n);
%             temp(5:n-1, 1:2) = cartesianCovarMatrix(6:n, 1:2);
%             temp(3:4, 5:n-1) = J * cartesianCovarMatrix(3:5, 6:n);
%             temp(5:n-1, 3:4) = cartesianCovarMatrix(6:n, 3:5) * J';
%             temp(5:n-1, 5:n-1) = cartesianCovarMatrix(6:n, 6:n);
%         end
% 
%         sphericalCovarMatrix = temp;
%     else
%         % If we have exactly a 5x5 matrix, output a 4x4 matrix
%         sphericalCovarMatrix = zeros(4, 4);
%         sphericalCovarMatrix(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
%         sphericalCovarMatrix(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
%         sphericalCovarMatrix(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
%         sphericalCovarMatrix(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
%     end
% end
% 
% function J = calculateJacobian(newangle1, newangle2, newangle3)
%     % Calculate the Jacobian matrix for transformation from
%     % Cartesian (newangle1, newangle2, newangle3) to Spherical (theta, phi)
%     %
%     % J = [∂θ/∂n1, ∂θ/∂n2, ∂θ/∂n3;
%     %      ∂φ/∂n1, ∂φ/∂n2, ∂φ/∂n3]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Partial derivatives for theta (acos(n3))
%     % ∂θ/∂n_i = -1/sqrt(1-n3^2) * ∂n3/∂n_i
% 
%     % Calculate ∂n3/∂n_i for each parameter
%     % For a normalized vector, ∂n_i/∂n_j involves the chain rule
%     den_theta = sqrt(1 - n3^2); % sqrt(1-cos²θ) = sin(θ)
% 
%     % Handle numerical issues near the poles
%     if abs(den_theta) < 1e-10
%         % Near poles (θ ≈ 0 or θ ≈ π), set derivatives to nominal values
%         dtheta_dn1 = 0;
%         dtheta_dn2 = 0;
% 
%         if n3 > 0  % Near north pole (θ ≈ 0)
%             dtheta_dn3 = -1;  % Approaching from positive side
%         else       % Near south pole (θ ≈ π)
%             dtheta_dn3 = 1;   % Approaching from negative side
%         end
%     else
%         % For components (n1,n2), we have:
%         % ∂n3/∂n_i = (∂/∂n_i)(n3/r) = (∂n3/∂n_i)/r - n3/(r^2) * (∂r/∂n_i)
%         % For n3 component itself:
%         % ∂n3/∂n3 = (∂/∂n3)(n3/r) = 1/r - n3/(r^2) * (n3/r)
% 
%         % Calculate ∂r/∂n_i for each parameter
%         dr_dn1 = n1; % ∂r/∂n1 = n1/r, but r has been normalized
%         dr_dn2 = n2; % ∂r/∂n2 = n2/r, but r has been normalized
%         dr_dn3 = n3; % ∂r/∂n3 = n3/r, but r has been normalized
% 
%         % Now calculate ∂n3/∂n_i
%         dn3_dn1 = -n3 * dr_dn1 / r;
%         dn3_dn2 = -n3 * dr_dn2 / r;
%         dn3_dn3 = 1/r - n3 * dr_dn3 / r;
% 
%         % Finally, calculate ∂θ/∂n_i
%         dtheta_dn1 = -dn3_dn1 / den_theta;
%         dtheta_dn2 = -dn3_dn2 / den_theta;
%         dtheta_dn3 = -dn3_dn3 / den_theta;
%     end
% 
%     % Partial derivatives for phi (atan2(n2, n1))
%     % ∂φ/∂n_i = (∂/∂n_i)atan2(n2, n1)
% 
%     % The derivative of atan2(y, x) is:
%     % ∂/∂x atan2(y, x) = -y/(x²+y²)
%     % ∂/∂y atan2(y, x) = x/(x²+y²)
% 
%     den_phi = n1^2 + n2^2;
% 
%     % Handle numerical issues when n1 and n2 are both near zero
%     if den_phi < 1e-10
%         % When close to the poles (n1≈0, n2≈0), phi becomes ill-defined
%         % Set derivatives to zero since phi changes have minimal effect near poles
%         dphi_dn1 = 0;
%         dphi_dn2 = 0;
%         dphi_dn3 = 0;
%     else
%         % Calculate normalized derivatives of phi
% 
%         % For components n1 and n2 directly involved in atan2:
%         % ∂φ/∂n1 = (∂/∂n1)atan2(n2, n1) = -n2/(n1²+n2²)
%         dphi_dn1_raw = -n2 / den_phi;
% 
%         % ∂φ/∂n2 = (∂/∂n2)atan2(n2, n1) = n1/(n1²+n2²)
%         dphi_dn2_raw = n1 / den_phi;
% 
%         % For n3, changes don't directly affect φ in the atan2 function,
%         % but they affect normalization:
%         dphi_dn3_raw = 0;
% 
%         % Now account for normalization effects on n1, n2
%         dn1_dn1 = 1/r - n1 * dr_dn1 / r;
%         dn1_dn2 = -n1 * dr_dn2 / r;
%         dn1_dn3 = -n1 * dr_dn3 / r;
% 
%         dn2_dn1 = -n2 * dr_dn1 / r;
%         dn2_dn2 = 1/r - n2 * dr_dn2 / r;
%         dn2_dn3 = -n2 * dr_dn3 / r;
% 
%         % Final derivatives using chain rule
%         dphi_dn1 = dphi_dn1_raw * dn1_dn1 + dphi_dn2_raw * dn2_dn1;
%         dphi_dn2 = dphi_dn1_raw * dn1_dn2 + dphi_dn2_raw * dn2_dn2;
%         dphi_dn3 = dphi_dn1_raw * dn1_dn3 + dphi_dn2_raw * dn2_dn3;
%     end
% 
%     % Construct the Jacobian matrix
%     J = [dtheta_dn1, dtheta_dn2, dtheta_dn3;
%          dphi_dn1, dphi_dn2, dphi_dn3];
% end
% 
% function currentFitPSF = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % This function calculates PSF using Cartesian dipole orientation representation
%     % Convert Cartesian to spherical representation for the models that require it
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % If it's a Gaussian model, we don't need orientations
%     if strcmpi(model, 'gaussian')
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
%     else
%         % For models that use dipole orientation
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
% 
%         % Convert Cartesian coordinates to spherical coordinates (theta, phi)
%         % assuming unit vector normalization
%         % Calculate theta and phi from Cartesian coordinates
%         normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%         % Normalize to ensure we have a unit vector
%         newangle1_norm = newangle1 / normFactor;
%         newangle2_norm = newangle2 / normFactor;
%         newangle3_norm = newangle3 / normFactor;
% 
%         % Calculate theta (inclination angle from z-axis)
%         theta = acos(newangle3_norm);
% 
%         % Calculate phi (azimuthal angle in xy-plane)
%         phi = atan2(newangle2_norm, newangle1_norm);
% 
%         % Ensure phi is in [0, 2π]
%         if phi < 0
%             phi = phi + 2*pi;
%         end
% 
%         paramEst.dipole = Dipole(theta, phi);
%     end
% 
%     if strcmpi(model, 'mortensen')
%         % Use Mortensen model
%         % Try two specific locations for the Python module (local vs cluster)
%         possibleDirs = {
%             fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
%             '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits' 
%         };
% 
%         % Try each location until we find the Python file
%         pyModuleFound = false;
%         for i = 1:length(possibleDirs)
%             pyDir = possibleDirs{i};
% 
%             % Check if directory exists
%             if ~exist(pyDir, 'dir')
%                 continue;  % Skip to next directory
%             end
% 
%             % Check if Python file exists in this directory
%             pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
%             if exist(pyFilePath, 'file')
%                 % Add to Python path
%                 if count(py.sys.path(), pyDir) == 0
%                     py.sys.path().insert(int32(0), pyDir);
%                 end
% 
%                 pyModuleFound = true;
%                 break;  % Found the file, stop looking
%             end
%         end
% 
%         % If the module wasn't found in any location, show an error
%         if ~pyModuleFound
%             error(['Python module "vectorized_mortensen_flipy_phi_theta.py" not found in any of the specified locations. ', ...
%                    'Please ensure the file exists in one of these directories:\n', ...
%                    '- %s\n', ...
%                    '- %s'], possibleDirs{1}, possibleDirs{2});
%         end
% 
%         % Run the function defined in your python file using the calculated theta and phi
%         currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
%             x, ... % x
%             y, ... % y
%             theta, ... % theta
%             phi, ... % phi
%             paramEst.nPixels,... % image_size_px
%             double(paramEst.pixelSize.inNanometer), ... % pixel_size_nm
%             double(paramEst.wavelength.inNanometer), ... % wavelength
%             paramEst.refractiveIndices(2), ... % n_objective
%             paramEst.refractiveIndices(1), ... % n_sample
%             paramEst.objectiveNA, ... % NA
%             photons ... % n_photons
%         );
% 
%         % Convert to matlab array
%         py_shape = currentPsf.shape;
%         rows = double(py_shape{1});
%         cols = double(py_shape{2});
% 
%         % Initialize MATLAB array
%         psf_matlab = zeros(rows, cols);
% 
%         % Copy values individually using item() method with explicit integer conversion
%         for i = 0:(rows-1)
%             for j = 0:(cols-1)
%                 % Use py.int to explicitly convert indices to Python integers
%                 psf_matlab(i+1, j+1) = double(currentPsf.item(py.tuple({py.int(i), py.int(j)})));
%             end
%         end
% 
%         currentPsf = psf_matlab;
% 
%         % this bit is the equivalent of the getIntensitiesCamera()
%         % stuff done in hinterer
%         totalIntensity = sum(sum(currentPsf));
%         currentPsf = currentPsf / totalIntensity * photons;
%         currentFitPSF = currentPsf;
% 
%     else
%         % Select appropriate BackFocalPlane function based on model
%         if strcmpi(model, 'gaussian')
%             bfp = BackFocalPlane_gaussian(paramEst); % Use Gaussian model
%         elseif strcmpi(model, 'hinterer')
%             bfp = BackFocalPlane(paramEst); % Use Hinterer model
%         end
% 
%         paramEst.backFocalPlane = bfp;
% 
%         % dave
%         paramEst.fieldBFP.x = bfp.electricField.x;
%         paramEst.fieldBFP.y = bfp.electricField.y;
% 
%         currentPsf = zeros(paramEst.nPixels,paramEst.nPixels); 
%         for k=1:size(paramEst.stageDrift.motion,1)
%             % Apply aberrations
%             aberrationCoeffs = getAberrations(paramEst,k);
%             fieldBFP = applyAberrations(paramEst, aberrationCoeffs);
%             % Get image from BFP field
%             currentPsf = currentPsf + getIntensitiesCamera(paramEst, fieldBFP);
%         end
% 
%         totalIntensity = sum(currentPsf,'all');
%         currentPsf = currentPsf ./ totalIntensity * photons + noise_est;
%         currentFitPSF = currentPsf;
%     end
% 
%     currentFitPSF = adjustExcitation(paramEst, currentFitPSF);
%     currentFitPSF = applyShotNoise(paramEst, currentFitPSF);
%     currentFitPSF = addBackgroundNoise(paramEst, currentFitPSF);
% end
% 
% function derivatives = calculateDerivativesCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % Calculate derivatives of the PSF with respect to each parameter including newangle1, newangle2, newangle3
%     % Using finite differences for approximation
% 
%     % Small step sizes for numerical derivatives
%     delta_pos = 1;     % nm
%     delta_angle = 0.01; % for Cartesian coordinates
% 
%     % Calculate PSF at current parameters
%     psf0 = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
% 
%     % Derivatives with respect to x and y position
%     psf_dx_plus = calculateProbabilityCartesian(paramEst, x + delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dx_minus = calculateProbabilityCartesian(paramEst, x - delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);
% 
%     psf_dy_plus = calculateProbabilityCartesian(paramEst, x, y + delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dy_minus = calculateProbabilityCartesian(paramEst, x, y - delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);
% 
%     % For models that only use position (like Gaussian)
%     if strcmpi(model, 'gaussian')
%         derivatives = {dx, dy};
%     else
%         % Derivative with respect to newangle1 (with boundary checks)
%         if newangle1 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf0 - psf_da1_minus) / delta_angle;
%         elseif newangle1 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf_da1_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle2 (with boundary checks)
%         if newangle2 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf0 - psf_da2_minus) / delta_angle;
%         elseif newangle2 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf_da2_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle3 (with boundary checks)
%         if newangle3 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf0 - psf_da3_minus) / delta_angle;
%         elseif newangle3 < delta_angle
%             % Near lower boundary, use forward difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf_da3_minus) / (2 * delta_angle);
%         end
% 
%         % Include orientation derivatives for dipole models
%         derivatives = {dx, dy, da1, da2, da3};
%     end
% end





% %% this version uses second-order derivs in jacobian, and works better
% 
% function [cartesianCovarMatrix, sphericalCovarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model)
%     % Calculates the theoretical covariance matrix for parameter estimates in Cartesian coordinates
%     % based on the Fisher Information Matrix, and optionally converts to spherical coordinates
%     %
%     % INPUTS:
%     %   paramEst - PSF object with image and parameters
%     %   x_est, y_est - Position estimates
%     %   defocus_est - Defocus estimate
%     %   newangle1_est, newangle2_est, newangle3_est - Cartesian dipole orientation estimates
%     %       newangle1 ∈ [-1,1], newangle2 ∈ [-1,1], newangle3 ∈ [0,1]
%     %   photon_est, noise_est - Photon and noise estimates
%     %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
%     %
%     % OUTPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (inverse of Fisher matrix)
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (theta, phi)
%     %   fisherMatrix - Fisher information matrix
% 
%     % Calculate probability distribution (PSF model)
%     psf = calculateProbabilityCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Calculate derivatives of PSF with respect to parameters
%     derivatives = calculateDerivativesCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Construct Fisher Information Matrix
%     fisherMatrix = zeros(length(derivatives), length(derivatives));
% 
%     % Fill Fisher Matrix
%     for i = 1:length(derivatives)
%         for j = 1:i  % Use symmetry to save calculations
%             % Formula for Fisher Matrix: I_ij = N * sum( (∂p/∂θ_i * ∂p/∂θ_j) / (p + b/N) )
%             % Where p is the probability, b is background, N is photon count
%             denom = psf + noise_est/photon_est;
%             fisherMatrix(i,j) = sum(sum(derivatives{i} .* derivatives{j} ./ denom));
% 
%             % Matrix is symmetric
%             if i ~= j
%                 fisherMatrix(j,i) = fisherMatrix(i,j);
%             end
%         end
%     end
% 
%     % Calculate Cartesian covariance matrix (Cramér-Rao lower bound)
%     cartesianCovarMatrix = inv(fisherMatrix);
% 
%     % Convert to spherical coordinates (if not using Gaussian model)
%     if ~strcmpi(model, 'gaussian') && nargout > 1
%         sphericalCovarMatrix = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1_est, newangle2_est, newangle3_est);
%     else
%         % For Gaussian model, just return the same covariance
%         sphericalCovarMatrix = cartesianCovarMatrix;
%     end
% end
% 
% function [sphericalCovarMatrix] = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1, newangle2, newangle3)
%     % Converts a covariance matrix from Cartesian coordinates (newangle1, newangle2, newangle3)
%     % to spherical coordinates (theta, phi) with second-order corrections
%     %
%     % INPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (5x5 or larger)
%     %                          Order is assumed to be [x, y, newangle1, newangle2, newangle3, ...]
%     %   newangle1, newangle2, newangle3 - Current values of Cartesian coordinates
%     %
%     % OUTPUTS:
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (4x4 or larger)
%     %                          Order is [x, y, theta, phi, ...]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize to ensure we have a unit vector
%     newangle1_norm = newangle1 / normFactor;
%     newangle2_norm = newangle2 / normFactor;
%     newangle3_norm = newangle3 / normFactor;
% 
%     % Calculate current theta and phi from the Cartesian coordinates
%     theta = acos(newangle3_norm);
%     phi = atan2(newangle2_norm, newangle1_norm);
% 
%     % Ensure phi is in [0, 2π]
%     if phi < 0
%         phi = phi + 2*pi;
%     end
% 
%     % Extract the orientation part of the covariance matrix (3x3 submatrix)
%     % This assumes the order in cartesianCovarMatrix is [x, y, newangle1, newangle2, newangle3, ...]
%     orientationCovar = cartesianCovarMatrix(3:5, 3:5);
% 
%     % Calculate the Jacobian matrix for the transformation (first-order terms)
%     % J_ij = ∂(spherical_i)/∂(cartesian_j)
%     J = calculateJacobian(newangle1, newangle2, newangle3);
% 
%     % Calculate second-order terms (Hessians)
%     H_theta = calculateThetaHessian(newangle1, newangle2, newangle3);
%     H_phi = calculatePhiHessian(newangle1, newangle2, newangle3);
% 
%     % First-order transformation of covariance matrix
%     % Cov_spherical = J * Cov_cartesian * J'
%     newOrientationCovar = J * orientationCovar * J';
% 
%     % Apply second-order corrections to the diagonal elements
%     theta_var_correction = 0.5 * trace(H_theta * orientationCovar);
%     phi_var_correction = 0.5 * trace(H_phi * orientationCovar);
% 
%     % Ensure corrections are reasonable (can be problematic near singularities)
%     if isnan(theta_var_correction) || isinf(theta_var_correction) || theta_var_correction < 0
%         theta_var_correction = 0;
%     end
% 
%     if isnan(phi_var_correction) || isinf(phi_var_correction) || phi_var_correction < 0
%         phi_var_correction = 0;
%     end
% 
%     % Apply corrections (only to the diagonal elements)
%     newOrientationCovar(1,1) = newOrientationCovar(1,1) + theta_var_correction;
%     newOrientationCovar(2,2) = newOrientationCovar(2,2) + phi_var_correction;
% 
%     % Create the new spherical covariance matrix, preserving other elements
%     sphericalCovarMatrix = cartesianCovarMatrix;
% 
%     % Replace the orientation part (3:5, 3:5) with the transformed 2x2 matrix (theta, phi)
%     % We'll need to adjust the matrix size since we're going from 3 parameters to 2
% 
%     % First, handle the case where we have a 5x5 or larger matrix
%     n = size(cartesianCovarMatrix, 1);
%     if n >= 5
%         % Create a new matrix of appropriate size (n-1 x n-1)
%         temp = zeros(n-1, n-1);
% 
%         % Copy x, y covariances directly
%         temp(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
% 
%         % Insert the transformed theta, phi covariances with second-order corrections
%         temp(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
% 
%         % Copy cross-covariances between x,y and orientation parameters
%         % For x,y vs theta,phi
%         temp(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
%         temp(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
% 
%         % Copy any remaining parameters (if matrix is larger than 5x5)
%         if n > 5
%             temp(1:2, 5:n-1) = cartesianCovarMatrix(1:2, 6:n);
%             temp(5:n-1, 1:2) = cartesianCovarMatrix(6:n, 1:2);
%             temp(3:4, 5:n-1) = J * cartesianCovarMatrix(3:5, 6:n);
%             temp(5:n-1, 3:4) = cartesianCovarMatrix(6:n, 3:5) * J';
%             temp(5:n-1, 5:n-1) = cartesianCovarMatrix(6:n, 6:n);
%         end
% 
%         sphericalCovarMatrix = temp;
%     else
%         % If we have exactly a 5x5 matrix, output a 4x4 matrix
%         sphericalCovarMatrix = zeros(4, 4);
%         sphericalCovarMatrix(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
%         sphericalCovarMatrix(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
%         sphericalCovarMatrix(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
%         sphericalCovarMatrix(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
%     end
% end
% 
% function H = calculateThetaHessian(newangle1, newangle2, newangle3)
%     % Calculate the Hessian matrix of theta with respect to (newangle1, newangle2, newangle3)
%     % H_ij = ∂²θ/∂n_i∂n_j
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);
%     newangle2 = max(min(newangle2, 1), -1);
%     newangle3 = max(min(newangle3, 1), 0);
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Intermediate terms for derivatives
%     den_theta = sqrt(1 - n3^2);
% 
%     % Handle numerical issues near the poles
%     if abs(den_theta) < 1e-10
%         % Near poles, return zero Hessian (approximation)
%         H = zeros(3, 3);
%         return;
%     end
% 
%     % Calculate partial derivatives of normalized components
%     dn1_dn1 = 1/r - n1^2/r;
%     dn1_dn2 = -n1*n2/r;
%     dn1_dn3 = -n1*n3/r;
% 
%     dn2_dn1 = -n2*n1/r;
%     dn2_dn2 = 1/r - n2^2/r;
%     dn2_dn3 = -n2*n3/r;
% 
%     dn3_dn1 = -n3*n1/r;
%     dn3_dn2 = -n3*n2/r;
%     dn3_dn3 = 1/r - n3^2/r;
% 
%     % First derivatives of theta
%     dtheta_dn3 = -1/den_theta;
% 
%     % Second derivatives - using chain rule and product rule
%     % ∂²θ/∂n_i∂n_j = ∂/∂n_i(∂θ/∂n_j)
% 
%     % For ∂²θ/∂n3² = ∂/∂n3(∂θ/∂n3)
%     d2theta_dn3_dn3 = -n3/(den_theta^3) * dn3_dn3;
% 
%     % For mixed derivatives ∂²θ/∂n_i∂n_3 = ∂/∂n_i(∂θ/∂n_3)
%     d2theta_dn1_dn3 = -n3/(den_theta^3) * dn3_dn1;
%     d2theta_dn2_dn3 = -n3/(den_theta^3) * dn3_dn2;
% 
%     % Other second derivatives are zero because ∂θ/∂n1 = ∂θ/∂n2 = 0 directly
%     d2theta_dn1_dn1 = 0;
%     d2theta_dn1_dn2 = 0;
%     d2theta_dn2_dn1 = 0;
%     d2theta_dn2_dn2 = 0;
% 
%     % Construct the Hessian matrix
%     H = [d2theta_dn1_dn1, d2theta_dn1_dn2, d2theta_dn1_dn3;
%          d2theta_dn2_dn1, d2theta_dn2_dn2, d2theta_dn2_dn3;
%          d2theta_dn1_dn3, d2theta_dn2_dn3, d2theta_dn3_dn3];
% end
% 
% function H = calculatePhiHessian(newangle1, newangle2, newangle3)
%     % Calculate the Hessian matrix of phi with respect to (newangle1, newangle2, newangle3)
%     % H_ij = ∂²φ/∂n_i∂n_j
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);
%     newangle2 = max(min(newangle2, 1), -1);
%     newangle3 = max(min(newangle3, 1), 0);
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Denominator for phi derivatives
%     den_phi = n1^2 + n2^2;
% 
%     % Handle numerical issues near the poles
%     if den_phi < 1e-10
%         % When close to the poles (n1≈0, n2≈0), phi derivatives become ill-defined
%         % Return zero Hessian as an approximation
%         H = zeros(3, 3);
%         return;
%     end
% 
%     % Calculate partial derivatives of normalized components
%     dn1_dn1 = 1/r - n1^2/r;
%     dn1_dn2 = -n1*n2/r;
%     dn1_dn3 = -n1*n3/r;
% 
%     dn2_dn1 = -n2*n1/r;
%     dn2_dn2 = 1/r - n2^2/r;
%     dn2_dn3 = -n2*n3/r;
% 
%     % First derivatives of phi
%     dphi_dn1 = -n2/den_phi;
%     dphi_dn2 = n1/den_phi;
% 
%     % Second derivatives using chain rule and product rule
% 
%     % ∂²φ/∂n1² = ∂/∂n1(∂φ/∂n1)
%     d2phi_dn1_dn1 = -dn2_dn1/den_phi + 2*n2*n1/den_phi^2 * dn1_dn1;
% 
%     % ∂²φ/∂n1∂n2 = ∂/∂n1(∂φ/∂n2)
%     d2phi_dn1_dn2 = dn1_dn2/den_phi - 2*n1*n2/den_phi^2 * dn1_dn2;
% 
%     % ∂²φ/∂n2∂n1 = ∂/∂n2(∂φ/∂n1)
%     d2phi_dn2_dn1 = -dn2_dn2/den_phi + 2*n2*n1/den_phi^2 * dn2_dn1;
% 
%     % ∂²φ/∂n2² = ∂/∂n2(∂φ/∂n2)
%     d2phi_dn2_dn2 = dn1_dn2/den_phi - 2*n1*n2/den_phi^2 * dn2_dn2;
% 
%     % ∂²φ/∂n3∂n1 and ∂²φ/∂n3∂n2
%     d2phi_dn3_dn1 = -dn2_dn3/den_phi + 2*n2*n1/den_phi^2 * dn1_dn3;
%     d2phi_dn3_dn2 = dn1_dn3/den_phi - 2*n1*n2/den_phi^2 * dn2_dn3;
% 
%     % ∂²φ/∂n3² = 0 because phi doesn't directly depend on n3
%     d2phi_dn3_dn3 = 0;
% 
%     % Construct the Hessian matrix
%     H = [d2phi_dn1_dn1, d2phi_dn1_dn2, d2phi_dn3_dn1;
%          d2phi_dn2_dn1, d2phi_dn2_dn2, d2phi_dn3_dn2;
%          d2phi_dn3_dn1, d2phi_dn3_dn2, d2phi_dn3_dn3];
% end
% 
% function J = calculateJacobian(newangle1, newangle2, newangle3)
%     % Calculate the Jacobian matrix for transformation from
%     % Cartesian (newangle1, newangle2, newangle3) to Spherical (theta, phi)
%     %
%     % J = [∂θ/∂n1, ∂θ/∂n2, ∂θ/∂n3;
%     %      ∂φ/∂n1, ∂φ/∂n2, ∂φ/∂n3]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Partial derivatives for theta (acos(n3))
%     % ∂θ/∂n_i = -1/sqrt(1-n3^2) * ∂n3/∂n_i
% 
%     % Calculate ∂n3/∂n_i for each parameter
%     % For a normalized vector, ∂n_i/∂n_j involves the chain rule
%     den_theta = sqrt(1 - n3^2); % sqrt(1-cos²θ) = sin(θ)
% 
%     % Handle numerical issues near the poles
%     if abs(den_theta) < 1e-10
%         % Near poles (θ ≈ 0 or θ ≈ π), set derivatives to nominal values
%         dtheta_dn1 = 0;
%         dtheta_dn2 = 0;
% 
%         if n3 > 0  % Near north pole (θ ≈ 0)
%             dtheta_dn3 = -1;  % Approaching from positive side
%         else       % Near south pole (θ ≈ π)
%             dtheta_dn3 = 1;   % Approaching from negative side
%         end
%     else
%         % For components (n1,n2), we have:
%         % ∂n3/∂n_i = (∂/∂n_i)(n3/r) = (∂n3/∂n_i)/r - n3/(r^2) * (∂r/∂n_i)
%         % For n3 component itself:
%         % ∂n3/∂n3 = (∂/∂n3)(n3/r) = 1/r - n3/(r^2) * (n3/r)
% 
%         % Calculate ∂r/∂n_i for each parameter
%         dr_dn1 = n1; % ∂r/∂n1 = n1/r, but r has been normalized
%         dr_dn2 = n2; % ∂r/∂n2 = n2/r, but r has been normalized
%         dr_dn3 = n3; % ∂r/∂n3 = n3/r, but r has been normalized
% 
%         % Now calculate ∂n3/∂n_i
%         dn3_dn1 = -n3 * dr_dn1 / r;
%         dn3_dn2 = -n3 * dr_dn2 / r;
%         dn3_dn3 = 1/r - n3 * dr_dn3 / r;
% 
%         % Finally, calculate ∂θ/∂n_i
%         dtheta_dn1 = -dn3_dn1 / den_theta;
%         dtheta_dn2 = -dn3_dn2 / den_theta;
%         dtheta_dn3 = -dn3_dn3 / den_theta;
%     end
% 
%     % Partial derivatives for phi (atan2(n2, n1))
%     % ∂φ/∂n_i = (∂/∂n_i)atan2(n2, n1)
% 
%     % The derivative of atan2(y, x) is:
%     % ∂/∂x atan2(y, x) = -y/(x²+y²)
%     % ∂/∂y atan2(y, x) = x/(x²+y²)
% 
%     den_phi = n1^2 + n2^2;
% 
%     % Handle numerical issues when n1 and n2 are both near zero
%     if den_phi < 1e-10
%         % When close to the poles (n1≈0, n2≈0), phi becomes ill-defined
%         % Set derivatives to zero since phi changes have minimal effect near poles
%         dphi_dn1 = 0;
%         dphi_dn2 = 0;
%         dphi_dn3 = 0;
%     else
%         % Calculate normalized derivatives of phi
% 
%         % For components n1 and n2 directly involved in atan2:
%         % ∂φ/∂n1 = (∂/∂n1)atan2(n2, n1) = -n2/(n1²+n2²)
%         dphi_dn1_raw = -n2 / den_phi;
% 
%         % ∂φ/∂n2 = (∂/∂n2)atan2(n2, n1) = n1/(n1²+n2²)
%         dphi_dn2_raw = n1 / den_phi;
% 
%         % For n3, changes don't directly affect φ in the atan2 function,
%         % but they affect normalization:
%         dphi_dn3_raw = 0;
% 
%         % Now account for normalization effects on n1, n2
%         dn1_dn1 = 1/r - n1 * dr_dn1 / r;
%         dn1_dn2 = -n1 * dr_dn2 / r;
%         dn1_dn3 = -n1 * dr_dn3 / r;
% 
%         dn2_dn1 = -n2 * dr_dn1 / r;
%         dn2_dn2 = 1/r - n2 * dr_dn2 / r;
%         dn2_dn3 = -n2 * dr_dn3 / r;
% 
%         % Final derivatives using chain rule
%         dphi_dn1 = dphi_dn1_raw * dn1_dn1 + dphi_dn2_raw * dn2_dn1;
%         dphi_dn2 = dphi_dn1_raw * dn1_dn2 + dphi_dn2_raw * dn2_dn2;
%         dphi_dn3 = dphi_dn1_raw * dn1_dn3 + dphi_dn2_raw * dn2_dn3;
%     end
% 
%     % Construct the Jacobian matrix
%     J = [dtheta_dn1, dtheta_dn2, dtheta_dn3;
%          dphi_dn1, dphi_dn2, dphi_dn3];
% end
% 
% function currentFitPSF = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % This function calculates PSF using Cartesian dipole orientation representation
%     % Convert Cartesian to spherical representation for the models that require it
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % If it's a Gaussian model, we don't need orientations
%     if strcmpi(model, 'gaussian')
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
%     else
%         % For models that use dipole orientation
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
% 
%         % Convert Cartesian coordinates to spherical coordinates (theta, phi)
%         % assuming unit vector normalization
%         % Calculate theta and phi from Cartesian coordinates
%         normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%         % Normalize to ensure we have a unit vector
%         newangle1_norm = newangle1 / normFactor;
%         newangle2_norm = newangle2 / normFactor;
%         newangle3_norm = newangle3 / normFactor;
% 
%         % Calculate theta (inclination angle from z-axis)
%         theta = acos(newangle3_norm);
% 
%         % Calculate phi (azimuthal angle in xy-plane)
%         phi = atan2(newangle2_norm, newangle1_norm);
% 
%         % Ensure phi is in [0, 2π]
%         if phi < 0
%             phi = phi + 2*pi;
%         end
% 
%         paramEst.dipole = Dipole(theta, phi);
%     end
% 
%     if strcmpi(model, 'mortensen')
%         % Use Mortensen model
%         % Try two specific locations for the Python module (local vs cluster)
%         possibleDirs = {
%             fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
%             '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits' 
%         };
% 
%         % Try each location until we find the Python file
%         pyModuleFound = false;
%         for i = 1:length(possibleDirs)
%             pyDir = possibleDirs{i};
% 
%             % Check if directory exists
%             if ~exist(pyDir, 'dir')
%                 continue;  % Skip to next directory
%             end
% 
%             % Check if Python file exists in this directory
%             pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
%             if exist(pyFilePath, 'file')
%                 % Add to Python path
%                 if count(py.sys.path(), pyDir) == 0
%                     py.sys.path().insert(int32(0), pyDir);
%                 end
% 
%                 pyModuleFound = true;
%                 break;  % Found the file, stop looking
%             end
%         end
% 
%         % If the module wasn't found in any location, show an error
%         if ~pyModuleFound
%             error(['Python module "vectorized_mortensen_flipy_phi_theta.py" not found in any of the specified locations. ', ...
%                    'Please ensure the file exists in one of these directories:\n', ...
%                    '- %s\n', ...
%                    '- %s'], possibleDirs{1}, possibleDirs{2});
%         end
% 
%         % Run the function defined in your python file using the calculated theta and phi
%         currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
%             x, ... % x
%             y, ... % y
%             theta, ... % theta
%             phi, ... % phi
%             paramEst.nPixels,... % image_size_px
%             double(paramEst.pixelSize.inNanometer), ... % pixel_size_nm
%             double(paramEst.wavelength.inNanometer), ... % wavelength
%             paramEst.refractiveIndices(2), ... % n_objective
%             paramEst.refractiveIndices(1), ... % n_sample
%             paramEst.objectiveNA, ... % NA
%             photons ... % n_photons
%         );
% 
%         % Convert to matlab array
%         py_shape = currentPsf.shape;
%         rows = double(py_shape{1});
%         cols = double(py_shape{2});
% 
%         % Initialize MATLAB array
%         psf_matlab = zeros(rows, cols);
% 
%         % Copy values individually using item() method with explicit integer conversion
%         for i = 0:(rows-1)
%             for j = 0:(cols-1)
%                 % Use py.int to explicitly convert indices to Python integers
%                 psf_matlab(i+1, j+1) = double(currentPsf.item(py.tuple({py.int(i), py.int(j)})));
%             end
%         end
% 
%         currentPsf = psf_matlab;
% 
%         % this bit is the equivalent of the getIntensitiesCamera()
%         % stuff done in hinterer
%         totalIntensity = sum(sum(currentPsf));
%         currentPsf = currentPsf / totalIntensity * photons;
%         currentFitPSF = currentPsf;
% 
%     else
%         % Select appropriate BackFocalPlane function based on model
%         if strcmpi(model, 'gaussian')
%             bfp = BackFocalPlane_gaussian(paramEst); % Use Gaussian model
%         elseif strcmpi(model, 'hinterer')
%             bfp = BackFocalPlane(paramEst); % Use Hinterer model
%         end
% 
%         paramEst.backFocalPlane = bfp;
% 
%         % dave
%         paramEst.fieldBFP.x = bfp.electricField.x;
%         paramEst.fieldBFP.y = bfp.electricField.y;
% 
%         currentPsf = zeros(paramEst.nPixels,paramEst.nPixels); 
%         for k=1:size(paramEst.stageDrift.motion,1)
%             % Apply aberrations
%             aberrationCoeffs = getAberrations(paramEst,k);
%             fieldBFP = applyAberrations(paramEst, aberrationCoeffs);
%             % Get image from BFP field
%             currentPsf = currentPsf + getIntensitiesCamera(paramEst, fieldBFP);
%         end
% 
%         totalIntensity = sum(currentPsf,'all');
%         currentPsf = currentPsf ./ totalIntensity * photons + noise_est;
%         currentFitPSF = currentPsf;
%     end
% 
%     currentFitPSF = adjustExcitation(paramEst, currentFitPSF);
%     currentFitPSF = applyShotNoise(paramEst, currentFitPSF);
%     currentFitPSF = addBackgroundNoise(paramEst, currentFitPSF);
% end
% 
% function derivatives = calculateDerivativesCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % Calculate derivatives of the PSF with respect to each parameter including newangle1, newangle2, newangle3
%     % Using finite differences for approximation
% 
%     % Small step sizes for numerical derivatives
%     delta_pos = 1;     % nm
%     delta_angle = 0.01; % for Cartesian coordinates
% 
%     % Calculate PSF at current parameters
%     psf0 = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
% 
%     % Derivatives with respect to x and y position
%     psf_dx_plus = calculateProbabilityCartesian(paramEst, x + delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dx_minus = calculateProbabilityCartesian(paramEst, x - delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);
% 
%     psf_dy_plus = calculateProbabilityCartesian(paramEst, x, y + delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dy_minus = calculateProbabilityCartesian(paramEst, x, y - delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);
% 
%     % For models that only use position (like Gaussian)
%     if strcmpi(model, 'gaussian')
%         derivatives = {dx, dy};
%     else
%         % Derivative with respect to newangle1 (with boundary checks)
%         if newangle1 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf0 - psf_da1_minus) / delta_angle;
%         elseif newangle1 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf_da1_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle2 (with boundary checks)
%         if newangle2 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf0 - psf_da2_minus) / delta_angle;
%         elseif newangle2 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf_da2_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle3 (with boundary checks)
%         if newangle3 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf0 - psf_da3_minus) / delta_angle;
%         elseif newangle3 < delta_angle
%             % Near lower boundary, use forward difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf_da3_minus) / (2 * delta_angle);
%         end
% 
%         % Include orientation derivatives for dipole models
%         derivatives = {dx, dy, da1, da2, da3};
%     end
% end




%% this version uses full riemannian thing as well as second-order stuff

function [cartesianCovarMatrix, sphericalCovarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model)
    % Calculates the theoretical covariance matrix for parameter estimates in Cartesian coordinates
    % based on the Fisher Information Matrix, and optionally converts to spherical coordinates
    %
    % INPUTS:
    %   paramEst - PSF object with image and parameters
    %   x_est, y_est - Position estimates
    %   defocus_est - Defocus estimate
    %   newangle1_est, newangle2_est, newangle3_est - Cartesian dipole orientation estimates
    %       newangle1 ∈ [-1,1], newangle2 ∈ [-1,1], newangle3 ∈ [0,1]
    %   photon_est, noise_est - Photon and noise estimates
    %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
    %
    % OUTPUTS:
    %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (inverse of Fisher matrix)
    %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (theta, phi)
    %   fisherMatrix - Fisher information matrix

    % Calculate probability distribution (PSF model)
    psf = calculateProbabilityCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);

    % Calculate derivatives of PSF with respect to parameters
    derivatives = calculateDerivativesCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);

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

    % % Calculate Cartesian covariance matrix (Cramér-Rao lower bound)
    % cartesianCovarMatrix = inv(fisherMatrix);

    % Set minimum eigenvalue to avoid singular inverse matrix - I tried
    % Tikhonov regularisation but this I guess is a better approach
    [V, D] = eig(fisherMatrix);
    eigenvalues = diag(D);
    max_eigenvalue = max(eigenvalues);
    
    % Set minimum eigenvalue threshold (e.g., 0.01% of maximum)
    min_eigenvalue = 1e-5 * max_eigenvalue;
    eigenvalues = max(eigenvalues, min_eigenvalue);
    
    % Reconstruct regularized fisher matrix
    reg_fisher = V * diag(eigenvalues) * V';
    cartesianCovarMatrix = inv(reg_fisher);



    % Convert to spherical coordinates (if not using Gaussian model)
    if ~strcmpi(model, 'gaussian') && nargout > 1
        sphericalCovarMatrix = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1_est, newangle2_est, newangle3_est);

        % Calculate the normalization factor
        normFactor = sqrt(newangle1_est^2 + newangle2_est^2 + newangle3_est^2);

        % Normalize to ensure we have a unit vector
        newangle3_norm = newangle3_est / normFactor;

        % Calculate theta (inclination angle from z-axis)
        theta = acos(newangle3_norm);

        % Apply Riemannian correction to spherical covariance matrix
        sphericalCovarMatrix = applyRiemannianCorrection(sphericalCovarMatrix, theta);
    else
        % For Gaussian model, just return the same covariance
        sphericalCovarMatrix = cartesianCovarMatrix;
    end
end

function [sphericalCovarMatrix] = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1, newangle2, newangle3)
    % Converts a covariance matrix from Cartesian coordinates (newangle1, newangle2, newangle3)
    % to spherical coordinates (theta, phi) with second-order corrections
    %
    % INPUTS:
    %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (5x5 or larger)
    %                          Order is assumed to be [x, y, newangle1, newangle2, newangle3, ...]
    %   newangle1, newangle2, newangle3 - Current values of Cartesian coordinates
    %
    % OUTPUTS:
    %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (4x4 or larger)
    %                          Order is [x, y, theta, phi, ...]

    % Ensure parameters are within their allowed ranges
    newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
    newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
    newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]

    % Calculate the normalization factor
    normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);

    % Normalize to ensure we have a unit vector
    newangle1_norm = newangle1 / normFactor;
    newangle2_norm = newangle2 / normFactor;
    newangle3_norm = newangle3 / normFactor;

    % Calculate current theta and phi from the Cartesian coordinates
    theta = acos(newangle3_norm);
    phi = atan2(newangle2_norm, newangle1_norm);

    % Ensure phi is in [0, 2π]
    if phi < 0
        phi = phi + 2*pi;
    end

    % Extract the orientation part of the covariance matrix (3x3 submatrix)
    % This assumes the order in cartesianCovarMatrix is [x, y, newangle1, newangle2, newangle3, ...]
    orientationCovar = cartesianCovarMatrix(3:5, 3:5);

    % Calculate the Jacobian matrix for the transformation (first-order terms)
    % J_ij = ∂(spherical_i)/∂(cartesian_j)
    J = calculateJacobian(newangle1, newangle2, newangle3);

    % Calculate second-order terms (Hessians)
    H_theta = calculateThetaHessian(newangle1, newangle2, newangle3);
    H_phi = calculatePhiHessian(newangle1, newangle2, newangle3);

    % First-order transformation of covariance matrix
    % Cov_spherical = J * Cov_cartesian * J'
    newOrientationCovar = J * orientationCovar * J';

    % Apply second-order corrections to the diagonal elements
    theta_var_correction = 0.5 * trace(H_theta * orientationCovar);
    phi_var_correction = 0.5 * trace(H_phi * orientationCovar);

    % Ensure corrections are reasonable (can be problematic near singularities)
    if isnan(theta_var_correction) || isinf(theta_var_correction) || theta_var_correction < 0
        theta_var_correction = 0;
    end

    if isnan(phi_var_correction) || isinf(phi_var_correction) || phi_var_correction < 0
        phi_var_correction = 0;
    end

    % Apply corrections (only to the diagonal elements)
    newOrientationCovar(1,1) = newOrientationCovar(1,1) + theta_var_correction;
    newOrientationCovar(2,2) = newOrientationCovar(2,2) + phi_var_correction;

    % Create the new spherical covariance matrix, preserving other elements
    sphericalCovarMatrix = cartesianCovarMatrix;

    % Replace the orientation part (3:5, 3:5) with the transformed 2x2 matrix (theta, phi)
    % We'll need to adjust the matrix size since we're going from 3 parameters to 2

    % First, handle the case where we have a 5x5 or larger matrix
    n = size(cartesianCovarMatrix, 1);
    if n >= 5
        % Create a new matrix of appropriate size (n-1 x n-1)
        temp = zeros(n-1, n-1);

        % Copy x, y covariances directly
        temp(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);

        % Insert the transformed theta, phi covariances with second-order corrections
        temp(3:4, 3:4) = newOrientationCovar(1:2, 1:2);

        % Copy cross-covariances between x,y and orientation parameters
        % For x,y vs theta,phi
        temp(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
        temp(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);

        % Copy any remaining parameters (if matrix is larger than 5x5)
        if n > 5
            temp(1:2, 5:n-1) = cartesianCovarMatrix(1:2, 6:n);
            temp(5:n-1, 1:2) = cartesianCovarMatrix(6:n, 1:2);
            temp(3:4, 5:n-1) = J * cartesianCovarMatrix(3:5, 6:n);
            temp(5:n-1, 3:4) = cartesianCovarMatrix(6:n, 3:5) * J';
            temp(5:n-1, 5:n-1) = cartesianCovarMatrix(6:n, 6:n);
        end

        sphericalCovarMatrix = temp;
    else
        % If we have exactly a 5x5 matrix, output a 4x4 matrix
        sphericalCovarMatrix = zeros(4, 4);
        sphericalCovarMatrix(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
        sphericalCovarMatrix(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
        sphericalCovarMatrix(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
        sphericalCovarMatrix(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
    end
end

function sphericalCovarMatrix = applyRiemannianCorrection(sphericalCovarMatrix, theta)
    % Apply Riemannian geometry correction to spherical covariance matrix
    % to account for the curved nature of the spherical parameter space
    %
    % INPUTS:
    %   sphericalCovarMatrix - Covariance matrix in spherical coordinates
    %   theta - Current value of theta (inclination angle)
    %
    % OUTPUTS:
    %   sphericalCovarMatrix - Corrected covariance matrix in spherical coordinates

    % Avoid singularities near the poles
    min_sin_theta = 1e-6;
    sin_theta = max(sin(theta), min_sin_theta);

    % Extract the angular part of the covariance matrix (theta, phi)
    angularCovar = sphericalCovarMatrix(3:4, 3:4);

    % The Riemannian metric tensor for spherical coordinates
    % g = [1, 0; 0, sin²(θ)]

    % Apply metric correction to the phi variance
    % This accounts for the fact that near poles, a small change in phi
    % corresponds to a much smaller physical change in direction than at the equator
    scaleFactor = 1 / (sin_theta^2);

    % Apply correction to phi variance (accounts for coordinate compression near poles)
    angularCovar(2,2) = angularCovar(2,2) * scaleFactor;

    % Also adjust the covariance between theta and phi
    angularCovar(1,2) = angularCovar(1,2) * sqrt(scaleFactor);
    angularCovar(2,1) = angularCovar(2,1) * sqrt(scaleFactor);

    % Update the angular part of the full covariance matrix
    sphericalCovarMatrix(3:4, 3:4) = angularCovar;

    % Adjust cross-covariances with other parameters
    n = size(sphericalCovarMatrix, 1);

    % Correct cross-covariances involving phi (row 4)
    for i = [1, 2, 5:n]
        if i <= n  % Make sure we're not exceeding matrix dimensions
            sphericalCovarMatrix(i,4) = sphericalCovarMatrix(i,4) * sqrt(scaleFactor);
            sphericalCovarMatrix(4,i) = sphericalCovarMatrix(4,i) * sqrt(scaleFactor);
        end
    end

    return;
end

function H = calculateThetaHessian(newangle1, newangle2, newangle3)
    % Calculate the Hessian matrix of theta with respect to (newangle1, newangle2, newangle3)
    % H_ij = ∂²θ/∂n_i∂n_j

    % Ensure parameters are within their allowed ranges
    newangle1 = max(min(newangle1, 1), -1);
    newangle2 = max(min(newangle2, 1), -1);
    newangle3 = max(min(newangle3, 1), 0);

    % Calculate the normalization factor
    r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);

    % Normalize components
    n1 = newangle1 / r;
    n2 = newangle2 / r;
    n3 = newangle3 / r;

    % Intermediate terms for derivatives
    den_theta = sqrt(1 - n3^2);

    % Handle numerical issues near the poles
    if abs(den_theta) < 1e-10
        % Near poles, return zero Hessian (approximation)
        H = zeros(3, 3);
        return;
    end

    % Calculate partial derivatives of normalized components
    dn1_dn1 = 1/r - n1^2/r;
    dn1_dn2 = -n1*n2/r;
    dn1_dn3 = -n1*n3/r;

    dn2_dn1 = -n2*n1/r;
    dn2_dn2 = 1/r - n2^2/r;
    dn2_dn3 = -n2*n3/r;

    dn3_dn1 = -n3*n1/r;
    dn3_dn2 = -n3*n2/r;
    dn3_dn3 = 1/r - n3^2/r;

    % First derivatives of theta
    dtheta_dn3 = -1/den_theta;

    % Second derivatives - using chain rule and product rule
    % ∂²θ/∂n_i∂n_j = ∂/∂n_i(∂θ/∂n_j)

    % For ∂²θ/∂n3² = ∂/∂n3(∂θ/∂n3)
    d2theta_dn3_dn3 = -n3/(den_theta^3) * dn3_dn3;

    % For mixed derivatives ∂²θ/∂n_i∂n_3 = ∂/∂n_i(∂θ/∂n_3)
    d2theta_dn1_dn3 = -n3/(den_theta^3) * dn3_dn1;
    d2theta_dn2_dn3 = -n3/(den_theta^3) * dn3_dn2;

    % Other second derivatives are zero because ∂θ/∂n1 = ∂θ/∂n2 = 0 directly
    d2theta_dn1_dn1 = 0;
    d2theta_dn1_dn2 = 0;
    d2theta_dn2_dn1 = 0;
    d2theta_dn2_dn2 = 0;

    % Construct the Hessian matrix
    H = [d2theta_dn1_dn1, d2theta_dn1_dn2, d2theta_dn1_dn3;
         d2theta_dn2_dn1, d2theta_dn2_dn2, d2theta_dn2_dn3;
         d2theta_dn1_dn3, d2theta_dn2_dn3, d2theta_dn3_dn3];
end

function H = calculatePhiHessian(newangle1, newangle2, newangle3)
    % Calculate the Hessian matrix of phi with respect to (newangle1, newangle2, newangle3)
    % H_ij = ∂²φ/∂n_i∂n_j

    % Ensure parameters are within their allowed ranges
    newangle1 = max(min(newangle1, 1), -1);
    newangle2 = max(min(newangle2, 1), -1);
    newangle3 = max(min(newangle3, 1), 0);

    % Calculate the normalization factor
    r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);

    % Normalize components
    n1 = newangle1 / r;
    n2 = newangle2 / r;
    n3 = newangle3 / r;

    % Denominator for phi derivatives
    den_phi = n1^2 + n2^2;

    % Handle numerical issues near the poles
    if den_phi < 1e-10
        % When close to the poles (n1≈0, n2≈0), phi derivatives become ill-defined
        % Return zero Hessian as an approximation
        H = zeros(3, 3);
        return;
    end

    % Calculate partial derivatives of normalized components
    dn1_dn1 = 1/r - n1^2/r;
    dn1_dn2 = -n1*n2/r;
    dn1_dn3 = -n1*n3/r;

    dn2_dn1 = -n2*n1/r;
    dn2_dn2 = 1/r - n2^2/r;
    dn2_dn3 = -n2*n3/r;

    % First derivatives of phi
    dphi_dn1 = -n2/den_phi;
    dphi_dn2 = n1/den_phi;

    % Second derivatives using chain rule and product rule

    % ∂²φ/∂n1² = ∂/∂n1(∂φ/∂n1)
    d2phi_dn1_dn1 = -dn2_dn1/den_phi + 2*n2*n1/den_phi^2 * dn1_dn1;

    % ∂²φ/∂n1∂n2 = ∂/∂n1(∂φ/∂n2)
    d2phi_dn1_dn2 = dn1_dn2/den_phi - 2*n1*n2/den_phi^2 * dn1_dn2;

    % ∂²φ/∂n2∂n1 = ∂/∂n2(∂φ/∂n1)
    d2phi_dn2_dn1 = -dn2_dn2/den_phi + 2*n2*n1/den_phi^2 * dn2_dn1;

    % ∂²φ/∂n2² = ∂/∂n2(∂φ/∂n2)
    d2phi_dn2_dn2 = dn1_dn2/den_phi - 2*n1*n2/den_phi^2 * dn2_dn2;

    % ∂²φ/∂n3∂n1 and ∂²φ/∂n3∂n2
    d2phi_dn3_dn1 = -dn2_dn3/den_phi + 2*n2*n1/den_phi^2 * dn1_dn3;
    d2phi_dn3_dn2 = dn1_dn3/den_phi - 2*n1*n2/den_phi^2 * dn2_dn3;

    % ∂²φ/∂n3² = 0 because phi doesn't directly depend on n3
    d2phi_dn3_dn3 = 0;

    % Construct the Hessian matrix
    H = [d2phi_dn1_dn1, d2phi_dn1_dn2, d2phi_dn3_dn1;
         d2phi_dn2_dn1, d2phi_dn2_dn2, d2phi_dn3_dn2;
         d2phi_dn3_dn1, d2phi_dn3_dn2, d2phi_dn3_dn3];
end

function J = calculateJacobian(newangle1, newangle2, newangle3)
    % Calculate the Jacobian matrix for transformation from
    % Cartesian (newangle1, newangle2, newangle3) to Spherical (theta, phi)
    %
    % J = [∂θ/∂n1, ∂θ/∂n2, ∂θ/∂n3;
    %      ∂φ/∂n1, ∂φ/∂n2, ∂φ/∂n3]

    % Ensure parameters are within their allowed ranges
    newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
    newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
    newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]

    % Calculate the normalization factor
    r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);

    % Normalize components
    n1 = newangle1 / r;
    n2 = newangle2 / r;
    n3 = newangle3 / r;

    % Partial derivatives for theta (acos(n3))
    % ∂θ/∂n_i = -1/sqrt(1-n3^2) * ∂n3/∂n_i

    % Calculate ∂n3/∂n_i for each parameter
    % For a normalized vector, ∂n_i/∂n_j involves the chain rule
    den_theta = sqrt(1 - n3^2); % sqrt(1-cos²θ) = sin(θ)

    % Handle numerical issues near the poles
    if abs(den_theta) < 1e-10
        % Near poles (θ ≈ 0 or θ ≈ π), set derivatives to nominal values
        dtheta_dn1 = 0;
        dtheta_dn2 = 0;

        if n3 > 0  % Near north pole (θ ≈ 0)
            dtheta_dn3 = -1;  % Approaching from positive side
        else       % Near south pole (θ ≈ π)
            dtheta_dn3 = 1;   % Approaching from negative side
        end
    else
        % For components (n1,n2), we have:
        % ∂n3/∂n_i = (∂/∂n_i)(n3/r) = (∂n3/∂n_i)/r - n3/(r^2) * (∂r/∂n_i)
        % For n3 component itself:
        % ∂n3/∂n3 = (∂/∂n3)(n3/r) = 1/r - n3/(r^2) * (n3/r)

        % Calculate ∂r/∂n_i for each parameter
        dr_dn1 = n1; % ∂r/∂n1 = n1/r, but r has been normalized
        dr_dn2 = n2; % ∂r/∂n2 = n2/r, but r has been normalized
        dr_dn3 = n3; % ∂r/∂n3 = n3/r, but r has been normalized

        % Now calculate ∂n3/∂n_i
        dn3_dn1 = -n3 * dr_dn1 / r;
        dn3_dn2 = -n3 * dr_dn2 / r;
        dn3_dn3 = 1/r - n3 * dr_dn3 / r;

        % Finally, calculate ∂θ/∂n_i
        dtheta_dn1 = -dn3_dn1 / den_theta;
        dtheta_dn2 = -dn3_dn2 / den_theta;
        dtheta_dn3 = -dn3_dn3 / den_theta;
    end

    % Partial derivatives for phi (atan2(n2, n1))
    % ∂φ/∂n_i = (∂/∂n_i)atan2(n2, n1)

    % The derivative of atan2(y, x) is:
    % ∂/∂x atan2(y, x) = -y/(x²+y²)
    % ∂/∂y atan2(y, x) = x/(x²+y²)

    den_phi = n1^2 + n2^2;

    % Handle numerical issues when n1 and n2 are both near zero
    if den_phi < 1e-10
        % When close to the poles (n1≈0, n2≈0), phi becomes ill-defined
        % Set derivatives to zero since phi changes have minimal effect near poles
        dphi_dn1 = 0;
        dphi_dn2 = 0;
        dphi_dn3 = 0;
    else
        % Calculate normalized derivatives of phi

        % For components n1 and n2 directly involved in atan2:
        % ∂φ/∂n1 = (∂/∂n1)atan2(n2, n1) = -n2/(n1²+n2²)
        dphi_dn1_raw = -n2 / den_phi;

        % ∂φ/∂n2 = (∂/∂n2)atan2(n2, n1) = n1/(n1²+n2²)
        dphi_dn2_raw = n1 / den_phi;

        % For n3, changes don't directly affect φ in the atan2 function,
        % but they affect normalization:
        dphi_dn3_raw = 0;

        % Now account for normalization effects on n1, n2
        dn1_dn1 = 1/r - n1 * dr_dn1 / r;
        dn1_dn2 = -n1 * dr_dn2 / r;
        dn1_dn3 = -n1 * dr_dn3 / r;

        dn2_dn1 = -n2 * dr_dn1 / r;
        dn2_dn2 = 1/r - n2 * dr_dn2 / r;
        dn2_dn3 = -n2 * dr_dn3 / r;

        % Final derivatives using chain rule
        dphi_dn1 = dphi_dn1_raw * dn1_dn1 + dphi_dn2_raw * dn2_dn1;
        dphi_dn2 = dphi_dn1_raw * dn1_dn2 + dphi_dn2_raw * dn2_dn2;
        dphi_dn3 = dphi_dn1_raw * dn1_dn3 + dphi_dn2_raw * dn2_dn3;
    end

    % Construct the Jacobian matrix
    J = [dtheta_dn1, dtheta_dn2, dtheta_dn3;
         dphi_dn1, dphi_dn2, dphi_dn3];
end

function currentFitPSF = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
    % This function calculates PSF using Cartesian dipole orientation representation
    % Convert Cartesian to spherical representation for the models that require it

    % Ensure parameters are within their allowed ranges
    newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
    newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
    newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]

    % If it's a Gaussian model, we don't need orientations
    if strcmpi(model, 'gaussian')
        paramEst.position = Length([x, y, 0], 'nm');
        paramEst.defocus = Length(defocus, 'nm');
    else
        % For models that use dipole orientation
        paramEst.position = Length([x, y, 0], 'nm');
        paramEst.defocus = Length(defocus, 'nm');

        % Convert Cartesian coordinates to spherical coordinates (theta, phi)
        % assuming unit vector normalization
        % Calculate theta and phi from Cartesian coordinates
        normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);

        % Normalize to ensure we have a unit vector
        newangle1_norm = newangle1 / normFactor;
        newangle2_norm = newangle2 / normFactor;
        newangle3_norm = newangle3 / normFactor;

        % Calculate theta (inclination angle from z-axis)
        theta = acos(newangle3_norm);

        % Calculate phi (azimuthal angle in xy-plane)
        phi = atan2(newangle2_norm, newangle1_norm);

        % Ensure phi is in [0, 2π]
        if phi < 0
            phi = phi + 2*pi;
        end

        paramEst.dipole = Dipole(theta, phi);
    end

    if strcmpi(model, 'mortensen')
        % Use Mortensen model
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

        % Run the function defined in your python file using the calculated theta and phi
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
        currentPsf = currentPsf ./ totalIntensity * photons + noise_est;
        currentFitPSF = currentPsf;
    end

    currentFitPSF = adjustExcitation(paramEst, currentFitPSF);
    currentFitPSF = applyShotNoise(paramEst, currentFitPSF);
    currentFitPSF = addBackgroundNoise(paramEst, currentFitPSF);
end

function derivatives = calculateDerivativesCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
    % Calculate derivatives of the PSF with respect to each parameter including newangle1, newangle2, newangle3
    % Using finite differences for approximation

    % Small step sizes for numerical derivatives
    delta_pos = 1;     % nm
    delta_angle = 0.01; % for Cartesian coordinates

    % Calculate PSF at current parameters
    psf0 = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);

    % Derivatives with respect to x and y position
    psf_dx_plus = calculateProbabilityCartesian(paramEst, x + delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    psf_dx_minus = calculateProbabilityCartesian(paramEst, x - delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);

    psf_dy_plus = calculateProbabilityCartesian(paramEst, x, y + delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    psf_dy_minus = calculateProbabilityCartesian(paramEst, x, y - delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
    dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);

    % For models that only use position (like Gaussian)
    if strcmpi(model, 'gaussian')
        derivatives = {dx, dy};
    else
        % Derivative with respect to newangle1 (with boundary checks)
        if newangle1 > 1 - delta_angle
            % Near upper boundary, use backward difference
            psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
            da1 = (psf0 - psf_da1_minus) / delta_angle;
        elseif newangle1 < -1 + delta_angle
            % Near lower boundary, use forward difference
            psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
            da1 = (psf_da1_plus - psf0) / delta_angle;
        else
            % Use central difference
            psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
            psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
            da1 = (psf_da1_plus - psf_da1_minus) / (2 * delta_angle);
        end

        % Derivative with respect to newangle2 (with boundary checks)
        if newangle2 > 1 - delta_angle
            % Near upper boundary, use backward difference
            psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
            da2 = (psf0 - psf_da2_minus) / delta_angle;
        elseif newangle2 < -1 + delta_angle
            % Near lower boundary, use forward difference
            psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
            da2 = (psf_da2_plus - psf0) / delta_angle;
        else
            % Use central difference
            psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
            psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
            da2 = (psf_da2_plus - psf_da2_minus) / (2 * delta_angle);
        end

        % Derivative with respect to newangle3 (with boundary checks)
        if newangle3 > 1 - delta_angle
            % Near upper boundary, use backward difference
            psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
            da3 = (psf0 - psf_da3_minus) / delta_angle;
        elseif newangle3 < delta_angle
            % Near lower boundary, use forward difference
            psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
            da3 = (psf_da3_plus - psf0) / delta_angle;
        else
            % Use central difference
            psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
            psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
            da3 = (psf_da3_plus - psf_da3_minus) / (2 * delta_angle);
        end

        % Include orientation derivatives for dipole models
        derivatives = {dx, dy, da1, da2, da3};
    end
end




% %% this version uses riemann, but not second-order
% 
% function [cartesianCovarMatrix, sphericalCovarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model)
%     % Calculates the theoretical covariance matrix for parameter estimates in Cartesian coordinates
%     % based on the Fisher Information Matrix, and optionally converts to spherical coordinates
%     %
%     % INPUTS:
%     %   paramEst - PSF object with image and parameters
%     %   x_est, y_est - Position estimates
%     %   defocus_est - Defocus estimate
%     %   newangle1_est, newangle2_est, newangle3_est - Cartesian dipole orientation estimates
%     %       newangle1 ∈ [-1,1], newangle2 ∈ [-1,1], newangle3 ∈ [0,1]
%     %   photon_est, noise_est - Photon and noise estimates
%     %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
%     %
%     % OUTPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (inverse of Fisher matrix)
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (theta, phi)
%     %   fisherMatrix - Fisher information matrix
% 
%     % Calculate probability distribution (PSF model)
%     psf = calculateProbabilityCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Calculate derivatives of PSF with respect to parameters
%     derivatives = calculateDerivativesCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Construct Fisher Information Matrix
%     fisherMatrix = zeros(length(derivatives), length(derivatives));
% 
%     % Fill Fisher Matrix
%     for i = 1:length(derivatives)
%         for j = 1:i  % Use symmetry to save calculations
%             % Formula for Fisher Matrix: I_ij = N * sum( (∂p/∂θ_i * ∂p/∂θ_j) / (p + b/N) )
%             % Where p is the probability, b is background, N is photon count
%             denom = psf + noise_est/photon_est;
%             fisherMatrix(i,j) = sum(sum(derivatives{i} .* derivatives{j} ./ denom));
% 
%             % Matrix is symmetric
%             if i ~= j
%                 fisherMatrix(j,i) = fisherMatrix(i,j);
%             end
%         end
%     end
% 
%     % Calculate Cartesian covariance matrix (Cramér-Rao lower bound)
%     cartesianCovarMatrix = inv(fisherMatrix);
% 
%     % Convert to spherical coordinates (if not using Gaussian model)
%     if ~strcmpi(model, 'gaussian') && nargout > 1
%         % Convert Cartesian to spherical covariance matrix
%         sphericalCovarMatrix = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1_est, newangle2_est, newangle3_est);
% 
%         % Calculate theta for Riemannian correction
%         normFactor = sqrt(newangle1_est^2 + newangle2_est^2 + newangle3_est^2);
%         newangle3_norm = newangle3_est / normFactor;
%         theta = acos(newangle3_norm);
% 
%         % Apply Riemannian correction to spherical covariance matrix
%         sphericalCovarMatrix = applyRiemannianCorrection(sphericalCovarMatrix, theta);
%     else
%         % For Gaussian model, just return the same covariance
%         sphericalCovarMatrix = cartesianCovarMatrix;
%     end
% end
% 
% function [sphericalCovarMatrix] = convertCartesianToSphericalCovariance(cartesianCovarMatrix, newangle1, newangle2, newangle3)
%     % Converts a covariance matrix from Cartesian coordinates (newangle1, newangle2, newangle3)
%     % to spherical coordinates (theta, phi)
%     %
%     % INPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (5x5 or larger)
%     %                          Order is assumed to be [x, y, newangle1, newangle2, newangle3, ...]
%     %   newangle1, newangle2, newangle3 - Current values of Cartesian coordinates
%     %
%     % OUTPUTS:
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (4x4 or larger)
%     %                          Order is [x, y, theta, phi, ...]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize to ensure we have a unit vector
%     newangle1_norm = newangle1 / normFactor;
%     newangle2_norm = newangle2 / normFactor;
%     newangle3_norm = newangle3 / normFactor;
% 
%     % Calculate current theta and phi from the Cartesian coordinates
%     theta = acos(newangle3_norm);
%     phi = atan2(newangle2_norm, newangle1_norm);
% 
%     % Ensure phi is in [0, 2π]
%     if phi < 0
%         phi = phi + 2*pi;
%     end
% 
%     % Extract the orientation part of the covariance matrix (3x3 submatrix)
%     % This assumes the order in cartesianCovarMatrix is [x, y, newangle1, newangle2, newangle3, ...]
%     orientationCovar = cartesianCovarMatrix(3:5, 3:5);
% 
%     % Calculate the Jacobian matrix for the transformation
%     % J_ij = ∂(spherical_i)/∂(cartesian_j)
%     J = calculateJacobian(newangle1, newangle2, newangle3);
% 
%     % Transform the orientation covariance using the Jacobian
%     % Cov_spherical = J * Cov_cartesian * J'
%     newOrientationCovar = J * orientationCovar * J';
% 
%     % Create the new spherical covariance matrix, preserving other elements
%     sphericalCovarMatrix = cartesianCovarMatrix;
% 
%     % Replace the orientation part (3:5, 3:5) with the transformed 2x2 matrix (theta, phi)
%     % We'll need to adjust the matrix size since we're going from 3 parameters to 2
% 
%     % First, handle the case where we have a 5x5 or larger matrix
%     n = size(cartesianCovarMatrix, 1);
%     if n >= 5
%         % Create a new matrix of appropriate size (n-1 x n-1)
%         temp = zeros(n-1, n-1);
% 
%         % Copy x, y covariances directly
%         temp(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
% 
%         % Insert the transformed theta, phi covariances
%         temp(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
% 
%         % Copy cross-covariances between x,y and orientation parameters
%         % For x,y vs theta,phi
%         temp(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
%         temp(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
% 
%         % Copy any remaining parameters (if matrix is larger than 5x5)
%         if n > 5
%             temp(1:2, 5:n-1) = cartesianCovarMatrix(1:2, 6:n);
%             temp(5:n-1, 1:2) = cartesianCovarMatrix(6:n, 1:2);
%             temp(3:4, 5:n-1) = J * cartesianCovarMatrix(3:5, 6:n);
%             temp(5:n-1, 3:4) = cartesianCovarMatrix(6:n, 3:5) * J';
%             temp(5:n-1, 5:n-1) = cartesianCovarMatrix(6:n, 6:n);
%         end
% 
%         sphericalCovarMatrix = temp;
%     else
%         % If we have exactly a 5x5 matrix, output a 4x4 matrix
%         sphericalCovarMatrix = zeros(4, 4);
%         sphericalCovarMatrix(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
%         sphericalCovarMatrix(3:4, 3:4) = newOrientationCovar(1:2, 1:2);
%         sphericalCovarMatrix(1:2, 3:4) = cartesianCovarMatrix(1:2, 3:5) * J';
%         sphericalCovarMatrix(3:4, 1:2) = J * cartesianCovarMatrix(3:5, 1:2);
%     end
% end
% 
% function sphericalCovarMatrix = applyRiemannianCorrection(sphericalCovarMatrix, theta)
%     % Apply Riemannian geometry correction to spherical covariance matrix
%     % to account for the curved nature of the spherical parameter space
%     %
%     % INPUTS:
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates
%     %   theta - Current value of theta (inclination angle)
%     %
%     % OUTPUTS:
%     %   sphericalCovarMatrix - Corrected covariance matrix in spherical coordinates
% 
%     % Avoid singularities near the poles
%     min_sin_theta = 1e-6;
%     sin_theta = max(sin(theta), min_sin_theta);
% 
%     % Extract the angular part of the covariance matrix (theta, phi)
%     angularCovar = sphericalCovarMatrix(3:4, 3:4);
% 
%     % The Riemannian metric tensor for spherical coordinates
%     % g = [1, 0; 0, sin²(θ)]
% 
%     % Apply metric correction to the phi variance
%     % This accounts for the fact that near poles, a small change in phi
%     % corresponds to a much smaller physical change in direction than at the equator
%     scaleFactor = 1 / (sin_theta^2);
% 
%     % Apply correction to phi variance (accounts for coordinate compression near poles)
%     angularCovar(2,2) = angularCovar(2,2) * scaleFactor;
% 
%     % Also adjust the covariance between theta and phi
%     angularCovar(1,2) = angularCovar(1,2) * sqrt(scaleFactor);
%     angularCovar(2,1) = angularCovar(2,1) * sqrt(scaleFactor);
% 
%     % Update the angular part of the full covariance matrix
%     sphericalCovarMatrix(3:4, 3:4) = angularCovar;
% 
%     % Adjust cross-covariances with other parameters
%     n = size(sphericalCovarMatrix, 1);
% 
%     % Correct cross-covariances involving phi (row 4)
%     for i = [1, 2, 5:n]
%         if i <= n  % Make sure we're not exceeding matrix dimensions
%             sphericalCovarMatrix(i,4) = sphericalCovarMatrix(i,4) * sqrt(scaleFactor);
%             sphericalCovarMatrix(4,i) = sphericalCovarMatrix(4,i) * sqrt(scaleFactor);
%         end
%     end
% end
% 
% function J = calculateJacobian(newangle1, newangle2, newangle3)
%     % Calculate the Jacobian matrix for transformation from
%     % Cartesian (newangle1, newangle2, newangle3) to Spherical (theta, phi)
%     %
%     % J = [∂θ/∂n1, ∂θ/∂n2, ∂θ/∂n3;
%     %      ∂φ/∂n1, ∂φ/∂n2, ∂φ/∂n3]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Partial derivatives for theta (acos(n3))
%     % ∂θ/∂n_i = -1/sqrt(1-n3^2) * ∂n3/∂n_i
% 
%     % Calculate ∂n3/∂n_i for each parameter
%     % For a normalized vector, ∂n_i/∂n_j involves the chain rule
%     den_theta = sqrt(1 - n3^2); % sqrt(1-cos²θ) = sin(θ)
% 
%     % Handle numerical issues near the poles
%     if abs(den_theta) < 1e-10
%         % Near poles (θ ≈ 0 or θ ≈ π), set derivatives to nominal values
%         dtheta_dn1 = 0;
%         dtheta_dn2 = 0;
% 
%         if n3 > 0  % Near north pole (θ ≈ 0)
%             dtheta_dn3 = -1;  % Approaching from positive side
%         else       % Near south pole (θ ≈ π)
%             dtheta_dn3 = 1;   % Approaching from negative side
%         end
%     else
%         % For components (n1,n2), we have:
%         % ∂n3/∂n_i = (∂/∂n_i)(n3/r) = (∂n3/∂n_i)/r - n3/(r^2) * (∂r/∂n_i)
%         % For n3 component itself:
%         % ∂n3/∂n3 = (∂/∂n3)(n3/r) = 1/r - n3/(r^2) * (n3/r)
% 
%         % Calculate ∂r/∂n_i for each parameter
%         dr_dn1 = n1; % ∂r/∂n1 = n1/r, but r has been normalized
%         dr_dn2 = n2; % ∂r/∂n2 = n2/r, but r has been normalized
%         dr_dn3 = n3; % ∂r/∂n3 = n3/r, but r has been normalized
% 
%         % Now calculate ∂n3/∂n_i
%         dn3_dn1 = -n3 * dr_dn1 / r;
%         dn3_dn2 = -n3 * dr_dn2 / r;
%         dn3_dn3 = 1/r - n3 * dr_dn3 / r;
% 
%         % Finally, calculate ∂θ/∂n_i
%         dtheta_dn1 = -dn3_dn1 / den_theta;
%         dtheta_dn2 = -dn3_dn2 / den_theta;
%         dtheta_dn3 = -dn3_dn3 / den_theta;
%     end
% 
%     % Partial derivatives for phi (atan2(n2, n1))
%     % ∂φ/∂n_i = (∂/∂n_i)atan2(n2, n1)
% 
%     % The derivative of atan2(y, x) is:
%     % ∂/∂x atan2(y, x) = -y/(x²+y²)
%     % ∂/∂y atan2(y, x) = x/(x²+y²)
% 
%     den_phi = n1^2 + n2^2;
% 
%     % Handle numerical issues when n1 and n2 are both near zero
%     if den_phi < 1e-10
%         % When close to the poles (n1≈0, n2≈0), phi becomes ill-defined
%         % Set derivatives to zero since phi changes have minimal effect near poles
%         dphi_dn1 = 0;
%         dphi_dn2 = 0;
%         dphi_dn3 = 0;
%     else
%         % Calculate normalized derivatives of phi
% 
%         % For components n1 and n2 directly involved in atan2:
%         % ∂φ/∂n1 = (∂/∂n1)atan2(n2, n1) = -n2/(n1²+n2²)
%         dphi_dn1_raw = -n2 / den_phi;
% 
%         % ∂φ/∂n2 = (∂/∂n2)atan2(n2, n1) = n1/(n1²+n2²)
%         dphi_dn2_raw = n1 / den_phi;
% 
%         % For n3, changes don't directly affect φ in the atan2 function,
%         % but they affect normalization:
%         dphi_dn3_raw = 0;
% 
%         % Now account for normalization effects on n1, n2
%         dn1_dn1 = 1/r - n1 * dr_dn1 / r;
%         dn1_dn2 = -n1 * dr_dn2 / r;
%         dn1_dn3 = -n1 * dr_dn3 / r;
% 
%         dn2_dn1 = -n2 * dr_dn1 / r;
%         dn2_dn2 = 1/r - n2 * dr_dn2 / r;
%         dn2_dn3 = -n2 * dr_dn3 / r;
% 
%         % Final derivatives using chain rule
%         dphi_dn1 = dphi_dn1_raw * dn1_dn1 + dphi_dn2_raw * dn2_dn1;
%         dphi_dn2 = dphi_dn1_raw * dn1_dn2 + dphi_dn2_raw * dn2_dn2;
%         dphi_dn3 = dphi_dn1_raw * dn1_dn3 + dphi_dn2_raw * dn2_dn3;
%     end
% 
%     % Construct the Jacobian matrix
%     J = [dtheta_dn1, dtheta_dn2, dtheta_dn3;
%          dphi_dn1, dphi_dn2, dphi_dn3];
% end
% 
% function currentFitPSF = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % This function calculates PSF using Cartesian dipole orientation representation
%     % Convert Cartesian to spherical representation for the models that require it
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % If it's a Gaussian model, we don't need orientations
%     if strcmpi(model, 'gaussian')
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
%     else
%         % For models that use dipole orientation
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
% 
%         % Convert Cartesian coordinates to spherical coordinates (theta, phi)
%         % assuming unit vector normalization
%         % Calculate theta and phi from Cartesian coordinates
%         normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%         % Normalize to ensure we have a unit vector
%         newangle1_norm = newangle1 / normFactor;
%         newangle2_norm = newangle2 / normFactor;
%         newangle3_norm = newangle3 / normFactor;
% 
%         % Calculate theta (inclination angle from z-axis)
%         theta = acos(newangle3_norm);
% 
%         % Calculate phi (azimuthal angle in xy-plane)
%         phi = atan2(newangle2_norm, newangle1_norm);
% 
%         % Ensure phi is in [0, 2π]
%         if phi < 0
%             phi = phi + 2*pi;
%         end
% 
%         paramEst.dipole = Dipole(theta, phi);
%     end
% 
%     if strcmpi(model, 'mortensen')
%         % Use Mortensen model
%         % Try two specific locations for the Python module (local vs cluster)
%         possibleDirs = {
%             fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
%             '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits' 
%         };
% 
%         % Try each location until we find the Python file
%         pyModuleFound = false;
%         for i = 1:length(possibleDirs)
%             pyDir = possibleDirs{i};
% 
%             % Check if directory exists
%             if ~exist(pyDir, 'dir')
%                 continue;  % Skip to next directory
%             end
% 
%             % Check if Python file exists in this directory
%             pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
%             if exist(pyFilePath, 'file')
%                 % Add to Python path
%                 if count(py.sys.path(), pyDir) == 0
%                     py.sys.path().insert(int32(0), pyDir);
%                 end
% 
%                 pyModuleFound = true;
%                 break;  % Found the file, stop looking
%             end
%         end
% 
%         % If the module wasn't found in any location, show an error
%         if ~pyModuleFound
%             error(['Python module "vectorized_mortensen_flipy_phi_theta.py" not found in any of the specified locations. ', ...
%                    'Please ensure the file exists in one of these directories:\n', ...
%                    '- %s\n', ...
%                    '- %s'], possibleDirs{1}, possibleDirs{2});
%         end
% 
%         % Run the function defined in your python file using the calculated theta and phi
%         currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
%             x, ... % x
%             y, ... % y
%             theta, ... % theta
%             phi, ... % phi
%             paramEst.nPixels,... % image_size_px
%             double(paramEst.pixelSize.inNanometer), ... % pixel_size_nm
%             double(paramEst.wavelength.inNanometer), ... % wavelength
%             paramEst.refractiveIndices(2), ... % n_objective
%             paramEst.refractiveIndices(1), ... % n_sample
%             paramEst.objectiveNA, ... % NA
%             photons ... % n_photons
%         );
% 
%         % Convert to matlab array
%         py_shape = currentPsf.shape;
%         rows = double(py_shape{1});
%         cols = double(py_shape{2});
% 
%         % Initialize MATLAB array
%         psf_matlab = zeros(rows, cols);
% 
%         % Copy values individually using item() method with explicit integer conversion
%         for i = 0:(rows-1)
%             for j = 0:(cols-1)
%                 % Use py.int to explicitly convert indices to Python integers
%                 psf_matlab(i+1, j+1) = double(currentPsf.item(py.tuple({py.int(i), py.int(j)})));
%             end
%         end
% 
%         currentPsf = psf_matlab;
% 
%         % this bit is the equivalent of the getIntensitiesCamera()
%         % stuff done in hinterer
%         totalIntensity = sum(sum(currentPsf));
%         currentPsf = currentPsf / totalIntensity * photons;
%         currentFitPSF = currentPsf;
% 
%     else
%         % Select appropriate BackFocalPlane function based on model
%         if strcmpi(model, 'gaussian')
%             bfp = BackFocalPlane_gaussian(paramEst); % Use Gaussian model
%         elseif strcmpi(model, 'hinterer')
%             bfp = BackFocalPlane(paramEst); % Use Hinterer model
%         end
% 
%         paramEst.backFocalPlane = bfp;
% 
%         % dave
%         paramEst.fieldBFP.x = bfp.electricField.x;
%         paramEst.fieldBFP.y = bfp.electricField.y;
% 
%         currentPsf = zeros(paramEst.nPixels,paramEst.nPixels); 
%         for k=1:size(paramEst.stageDrift.motion,1)
%             % Apply aberrations
%             aberrationCoeffs = getAberrations(paramEst,k);
%             fieldBFP = applyAberrations(paramEst, aberrationCoeffs);
%             % Get image from BFP field
%             currentPsf = currentPsf + getIntensitiesCamera(paramEst, fieldBFP);
%         end
% 
%         totalIntensity = sum(currentPsf,'all');
%         currentPsf = currentPsf ./ totalIntensity * photons + noise_est;
%         currentFitPSF = currentPsf;
%     end
% 
%     currentFitPSF = adjustExcitation(paramEst, currentFitPSF);
%     currentFitPSF = applyShotNoise(paramEst, currentFitPSF);
%     currentFitPSF = addBackgroundNoise(paramEst, currentFitPSF);
% end
% 
% function derivatives = calculateDerivativesCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % Calculate derivatives of the PSF with respect to each parameter including newangle1, newangle2, newangle3
%     % Using finite differences for approximation
% 
%     % Small step sizes for numerical derivatives
%     delta_pos = 1;     % nm
%     delta_angle = 0.01; % for Cartesian coordinates
% 
%     % Calculate PSF at current parameters
%     psf0 = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
% 
%     % Derivatives with respect to x and y position
%     psf_dx_plus = calculateProbabilityCartesian(paramEst, x + delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dx_minus = calculateProbabilityCartesian(paramEst, x - delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);
% 
%     psf_dy_plus = calculateProbabilityCartesian(paramEst, x, y + delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dy_minus = calculateProbabilityCartesian(paramEst, x, y - delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);
% 
%     % For models that only use position (like Gaussian)
%     if strcmpi(model, 'gaussian')
%         derivatives = {dx, dy};
%     else
%         % Derivative with respect to newangle1 (with boundary checks)
%         if newangle1 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf0 - psf_da1_minus) / delta_angle;
%         elseif newangle1 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf_da1_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle2 (with boundary checks)
%         if newangle2 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf0 - psf_da2_minus) / delta_angle;
%         elseif newangle2 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf_da2_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle3 (with boundary checks)
%         if newangle3 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf0 - psf_da3_minus) / delta_angle;
%         elseif newangle3 < delta_angle
%             % Near lower boundary, use forward difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf_da3_minus) / (2 * delta_angle);
%         end
% 
%         % Include orientation derivatives for dipole models
%         derivatives = {dx, dy, da1, da2, da3};
%     end
% end

% %% this uses riemann and second-order, but I think without the redundancy that the earlier version had
% 
% function [cartesianCovarMatrix, sphericalCovarMatrix, fisherMatrix] = calculateCartesianCoordinatesCovarianceMatrix(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model)
%     % Calculates the theoretical covariance matrix for parameter estimates in Cartesian coordinates
%     % based on the Fisher Information Matrix, and optionally converts to spherical coordinates
%     %
%     % INPUTS:
%     %   paramEst - PSF object with image and parameters
%     %   x_est, y_est - Position estimates
%     %   defocus_est - Defocus estimate
%     %   newangle1_est, newangle2_est, newangle3_est - Cartesian dipole orientation estimates
%     %       newangle1 ∈ [-1,1], newangle2 ∈ [-1,1], newangle3 ∈ [0,1]
%     %   photon_est, noise_est - Photon and noise estimates
%     %   model - PSF model being used ('mortensen', 'hinterer', or 'gaussian')
%     %
%     % OUTPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (inverse of Fisher matrix)
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (theta, phi)
%     %   fisherMatrix - Fisher information matrix
% 
%     % Calculate probability distribution (PSF model)
%     psf = calculateProbabilityCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Calculate derivatives of PSF with respect to parameters
%     derivatives = calculateDerivativesCartesian(paramEst, x_est, y_est, defocus_est, newangle1_est, newangle2_est, newangle3_est, photon_est, noise_est, model);
% 
%     % Construct Fisher Information Matrix
%     fisherMatrix = zeros(length(derivatives), length(derivatives));
% 
%     % Fill Fisher Matrix
%     for i = 1:length(derivatives)
%         for j = 1:i  % Use symmetry to save calculations
%             % Formula for Fisher Matrix: I_ij = N * sum( (∂p/∂θ_i * ∂p/∂θ_j) / (p + b/N) )
%             % Where p is the probability, b is background, N is photon count
%             denom = psf + noise_est/photon_est;
%             fisherMatrix(i,j) = sum(sum(derivatives{i} .* derivatives{j} ./ denom));
% 
%             % Matrix is symmetric
%             if i ~= j
%                 fisherMatrix(j,i) = fisherMatrix(i,j);
%             end
%         end
%     end
% 
%     % Calculate Cartesian covariance matrix (Cramér-Rao lower bound)
%     cartesianCovarMatrix = inv(fisherMatrix);
% 
%     % Convert to spherical coordinates (if not using Gaussian model)
%     if ~strcmpi(model, 'gaussian') && nargout > 1
%         % Calculate normalization and theta for later use
%         normFactor = sqrt(newangle1_est^2 + newangle2_est^2 + newangle3_est^2);
%         newangle3_norm = newangle3_est / normFactor;
%         theta = acos(newangle3_norm);
% 
%         % Convert from Cartesian to spherical coordinates using theoretically sound approach
%         sphericalCovarMatrix = convertCartesianToSphericalCovarianceComplete( ...
%             cartesianCovarMatrix, newangle1_est, newangle2_est, newangle3_est, theta);
%     else
%         % For Gaussian model, just return the same covariance
%         sphericalCovarMatrix = cartesianCovarMatrix;
%     end
% end
% 
% function [sphericalCovarMatrix] = convertCartesianToSphericalCovarianceComplete( ...
%     cartesianCovarMatrix, newangle1, newangle2, newangle3, theta)
%     % Converts a covariance matrix from Cartesian coordinates to spherical coordinates
%     % with a complete theoretical approach incorporating both second-order (Hessian) terms
%     % and Riemannian metric corrections. This follows the full transformation law for
%     % probability distributions under nonlinear coordinate transformations.
%     %
%     % INPUTS:
%     %   cartesianCovarMatrix - Covariance matrix in Cartesian coordinates (5x5 or larger)
%     %                          Order is assumed to be [x, y, newangle1, newangle2, newangle3, ...]
%     %   newangle1, newangle2, newangle3 - Current values of Cartesian coordinates
%     %   theta - Pre-calculated theta value
%     %
%     % OUTPUTS:
%     %   sphericalCovarMatrix - Covariance matrix in spherical coordinates (4x4 or larger)
%     %                          Order is [x, y, theta, phi, ...]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);
%     newangle2 = max(min(newangle2, 1), -1);
%     newangle3 = max(min(newangle3, 1), 0);
% 
%     % Calculate the normalization factor
%     normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize to ensure we have a unit vector
%     newangle1_norm = newangle1 / normFactor;
%     newangle2_norm = newangle2 / normFactor;
%     newangle3_norm = newangle3 / normFactor;
% 
%     % Calculate phi from the Cartesian coordinates
%     phi = atan2(newangle2_norm, newangle1_norm);
%     if phi < 0
%         phi = phi + 2*pi;
%     end
% 
%     % Extract the orientation part of the covariance matrix (3x3 submatrix)
%     orientationCovar = cartesianCovarMatrix(3:5, 3:5);
% 
%     %------------------------------------------------------------------
%     % Step 1: First-order transformation using the Jacobian
%     %------------------------------------------------------------------
% 
%     % Calculate the Jacobian matrix for the transformation
%     J = calculateJacobian(newangle1, newangle2, newangle3);
% 
%     % First-order transformation of covariance matrix
%     % Σ_spherical_first_order = J * Σ_cartesian * J'
%     sphericalOrientationCovar_1st = J * orientationCovar * J';
% 
%     %------------------------------------------------------------------
%     % Step 2: Add second-order corrections (Hessian terms)
%     %------------------------------------------------------------------
% 
%     % Calculate Hessian matrices for theta and phi
%     H_theta = calculateThetaHessian(newangle1, newangle2, newangle3);
%     H_phi = calculatePhiHessian(newangle1, newangle2, newangle3);
% 
%     % Calculate second-order corrections for variance terms
%     % These account for the nonlinearity of the coordinate transformation
%     % The formula is: correction = 0.5 * trace(H * Σ)
%     theta_var_correction = 0.5 * trace(H_theta * orientationCovar);
%     phi_var_correction = 0.5 * trace(H_phi * orientationCovar);
% 
%     % Ensure corrections are numerically valid
%     if isnan(theta_var_correction) || isinf(theta_var_correction) || theta_var_correction < 0
%         theta_var_correction = 0;
%     end
% 
%     if isnan(phi_var_correction) || isinf(phi_var_correction) || phi_var_correction < 0
%         phi_var_correction = 0;
%     end
% 
%     % Apply second-order corrections to the diagonal elements of the covariance
%     sphericalOrientationCovar_2nd = sphericalOrientationCovar_1st;
%     sphericalOrientationCovar_2nd(1,1) = sphericalOrientationCovar_1st(1,1) + theta_var_correction;
%     sphericalOrientationCovar_2nd(2,2) = sphericalOrientationCovar_1st(2,2) + phi_var_correction;
% 
%     %------------------------------------------------------------------
%     % Step 3: Apply Riemannian metric correction
%     %------------------------------------------------------------------
% 
%     % Avoid singularities near the poles
%     min_sin_theta = 1e-6;
%     sin_theta = max(sin(theta), min_sin_theta);
% 
%     % Create the inverse metric tensor G^(-1)
%     % For spherical coordinates on a unit sphere, G^(-1) = [1, 0; 0, 1/sin(θ)]
%     G_inv = [1, 0; 0, 1/sin_theta];
% 
%     % Apply the Riemannian correction
%     % Σ_riemannian = G^(-1) * Σ_spherical * G^(-1)
%     riemannianOrientationCovar = G_inv * sphericalOrientationCovar_2nd * G_inv;
% 
%     %------------------------------------------------------------------
%     % Step 4: Create the complete transformed covariance matrix
%     %------------------------------------------------------------------
% 
%     % First, handle the case where we have a 5x5 or larger matrix
%     n = size(cartesianCovarMatrix, 1);
%     if n >= 5
%         % Create a new matrix of appropriate size (n-1 x n-1)
%         temp = zeros(n-1, n-1);
% 
%         % Copy x, y covariances directly
%         temp(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
% 
%         % Insert the fully transformed theta, phi covariances
%         temp(3:4, 3:4) = riemannianOrientationCovar;
% 
%         % Transform cross-covariances between x,y and orientation parameters
%         % First, apply Jacobian transformation
%         xy_orientation_cross = cartesianCovarMatrix(1:2, 3:5) * J';
% 
%         % Then apply the Riemannian metric to the phi component
%         xy_orientation_cross(:, 2) = xy_orientation_cross(:, 2) / sin_theta;
% 
%         temp(1:2, 3:4) = xy_orientation_cross;
% 
%         % Similarly for transpose cross-covariances
%         orientation_xy_cross = J * cartesianCovarMatrix(3:5, 1:2);
%         orientation_xy_cross(2, :) = orientation_xy_cross(2, :) / sin_theta;
% 
%         temp(3:4, 1:2) = orientation_xy_cross;
% 
%         % Copy any remaining parameters (if matrix is larger than 5x5)
%         if n > 5
%             temp(1:2, 5:n-1) = cartesianCovarMatrix(1:2, 6:n);
%             temp(5:n-1, 1:2) = cartesianCovarMatrix(6:n, 1:2);
% 
%             % Transform cross-covariances involving orientation parameters
%             orientation_other_cross = J * cartesianCovarMatrix(3:5, 6:n);
%             orientation_other_cross(2, :) = orientation_other_cross(2, :) / sin_theta;
%             temp(3:4, 5:n-1) = orientation_other_cross;
% 
%             other_orientation_cross = cartesianCovarMatrix(6:n, 3:5) * J';
%             other_orientation_cross(:, 2) = other_orientation_cross(:, 2) / sin_theta;
%             temp(5:n-1, 3:4) = other_orientation_cross;
% 
%             temp(5:n-1, 5:n-1) = cartesianCovarMatrix(6:n, 6:n);
%         end
% 
%         sphericalCovarMatrix = temp;
%     else
%         % If we have exactly a 5x5 matrix, output a 4x4 matrix
%         sphericalCovarMatrix = zeros(4, 4);
%         sphericalCovarMatrix(1:2, 1:2) = cartesianCovarMatrix(1:2, 1:2);
%         sphericalCovarMatrix(3:4, 3:4) = riemannianOrientationCovar;
% 
%         % Transform cross-covariances
%         xy_orientation_cross = cartesianCovarMatrix(1:2, 3:5) * J';
%         xy_orientation_cross(:, 2) = xy_orientation_cross(:, 2) / sin_theta;
%         sphericalCovarMatrix(1:2, 3:4) = xy_orientation_cross;
% 
%         orientation_xy_cross = J * cartesianCovarMatrix(3:5, 1:2);
%         orientation_xy_cross(2, :) = orientation_xy_cross(2, :) / sin_theta;
%         sphericalCovarMatrix(3:4, 1:2) = orientation_xy_cross;
%     end
% end
% 
% function H = calculateThetaHessian(newangle1, newangle2, newangle3)
%     % Calculate the Hessian matrix of theta with respect to (newangle1, newangle2, newangle3)
%     % H_ij = ∂²θ/∂n_i∂n_j
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);
%     newangle2 = max(min(newangle2, 1), -1);
%     newangle3 = max(min(newangle3, 1), 0);
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Intermediate terms for derivatives
%     den_theta = sqrt(1 - n3^2);
% 
%     % Handle numerical issues near the poles
%     if abs(den_theta) < 1e-10
%         % Near poles, return zero Hessian (approximation)
%         H = zeros(3, 3);
%         return;
%     end
% 
%     % Calculate partial derivatives of normalized components
%     dn1_dn1 = 1/r - n1^2/r;
%     dn1_dn2 = -n1*n2/r;
%     dn1_dn3 = -n1*n3/r;
% 
%     dn2_dn1 = -n2*n1/r;
%     dn2_dn2 = 1/r - n2^2/r;
%     dn2_dn3 = -n2*n3/r;
% 
%     dn3_dn1 = -n3*n1/r;
%     dn3_dn2 = -n3*n2/r;
%     dn3_dn3 = 1/r - n3^2/r;
% 
%     % First derivatives of theta
%     dtheta_dn3 = -1/den_theta;
% 
%     % Second derivatives - using chain rule and product rule
%     % ∂²θ/∂n_i∂n_j = ∂/∂n_i(∂θ/∂n_j)
% 
%     % For ∂²θ/∂n3² = ∂/∂n3(∂θ/∂n3)
%     d2theta_dn3_dn3 = -n3/(den_theta^3) * dn3_dn3;
% 
%     % For mixed derivatives ∂²θ/∂n_i∂n_3 = ∂/∂n_i(∂θ/∂n_3)
%     d2theta_dn1_dn3 = -n3/(den_theta^3) * dn3_dn1;
%     d2theta_dn2_dn3 = -n3/(den_theta^3) * dn3_dn2;
% 
%     % Other second derivatives are zero because ∂θ/∂n1 = ∂θ/∂n2 = 0 directly
%     d2theta_dn1_dn1 = 0;
%     d2theta_dn1_dn2 = 0;
%     d2theta_dn2_dn1 = 0;
%     d2theta_dn2_dn2 = 0;
% 
%     % Construct the Hessian matrix
%     H = [d2theta_dn1_dn1, d2theta_dn1_dn2, d2theta_dn1_dn3;
%          d2theta_dn2_dn1, d2theta_dn2_dn2, d2theta_dn2_dn3;
%          d2theta_dn1_dn3, d2theta_dn2_dn3, d2theta_dn3_dn3];
% end
% 
% function H = calculatePhiHessian(newangle1, newangle2, newangle3)
%     % Calculate the Hessian matrix of phi with respect to (newangle1, newangle2, newangle3)
%     % H_ij = ∂²φ/∂n_i∂n_j
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);
%     newangle2 = max(min(newangle2, 1), -1);
%     newangle3 = max(min(newangle3, 1), 0);
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Denominator for phi derivatives
%     den_phi = n1^2 + n2^2;
% 
%     % Handle numerical issues near the poles
%     if den_phi < 1e-10
%         % When close to the poles (n1≈0, n2≈0), phi derivatives become ill-defined
%         % Return zero Hessian as an approximation
%         H = zeros(3, 3);
%         return;
%     end
% 
%     % Calculate partial derivatives of normalized components
%     dn1_dn1 = 1/r - n1^2/r;
%     dn1_dn2 = -n1*n2/r;
%     dn1_dn3 = -n1*n3/r;
% 
%     dn2_dn1 = -n2*n1/r;
%     dn2_dn2 = 1/r - n2^2/r;
%     dn2_dn3 = -n2*n3/r;
% 
%     % First derivatives of phi
%     dphi_dn1 = -n2/den_phi;
%     dphi_dn2 = n1/den_phi;
% 
%     % Second derivatives using chain rule and product rule
% 
%     % ∂²φ/∂n1² = ∂/∂n1(∂φ/∂n1)
%     d2phi_dn1_dn1 = -dn2_dn1/den_phi + 2*n2*n1/den_phi^2 * dn1_dn1;
% 
%     % ∂²φ/∂n1∂n2 = ∂/∂n1(∂φ/∂n2)
%     d2phi_dn1_dn2 = dn1_dn2/den_phi - 2*n1*n2/den_phi^2 * dn1_dn2;
% 
%     % ∂²φ/∂n2∂n1 = ∂/∂n2(∂φ/∂n1)
%     d2phi_dn2_dn1 = -dn2_dn2/den_phi + 2*n2*n1/den_phi^2 * dn2_dn1;
% 
%     % ∂²φ/∂n2² = ∂/∂n2(∂φ/∂n2)
%     d2phi_dn2_dn2 = dn1_dn2/den_phi - 2*n1*n2/den_phi^2 * dn2_dn2;
% 
%     % ∂²φ/∂n3∂n1 and ∂²φ/∂n3∂n2
%     d2phi_dn3_dn1 = -dn2_dn3/den_phi + 2*n2*n1/den_phi^2 * dn1_dn3;
%     d2phi_dn3_dn2 = dn1_dn3/den_phi - 2*n1*n2/den_phi^2 * dn2_dn3;
% 
%     % ∂²φ/∂n3² = 0 because phi doesn't directly depend on n3
%     d2phi_dn3_dn3 = 0;
% 
%     % Construct the Hessian matrix
%     H = [d2phi_dn1_dn1, d2phi_dn1_dn2, d2phi_dn3_dn1;
%          d2phi_dn2_dn1, d2phi_dn2_dn2, d2phi_dn3_dn2;
%          d2phi_dn3_dn1, d2phi_dn3_dn2, d2phi_dn3_dn3];
% end
% 
% function J = calculateJacobian(newangle1, newangle2, newangle3)
%     % Calculate the Jacobian matrix for transformation from
%     % Cartesian (newangle1, newangle2, newangle3) to Spherical (theta, phi)
%     %
%     % J = [∂θ/∂n1, ∂θ/∂n2, ∂θ/∂n3;
%     %      ∂φ/∂n1, ∂φ/∂n2, ∂φ/∂n3]
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % Calculate the normalization factor
%     r = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%     % Normalize components
%     n1 = newangle1 / r;
%     n2 = newangle2 / r;
%     n3 = newangle3 / r;
% 
%     % Partial derivatives for theta (acos(n3))
%     % ∂θ/∂n_i = -1/sqrt(1-n3^2) * ∂n3/∂n_i
% 
%     % Calculate ∂n3/∂n_i for each parameter
%     % For a normalized vector, ∂n_i/∂n_j involves the chain rule
%     den_theta = sqrt(1 - n3^2); % sqrt(1-cos²θ) = sin(θ)
% 
%     % Handle numerical issues near the poles
%     if abs(den_theta) < 1e-10
%         % Near poles (θ ≈ 0 or θ ≈ π), set derivatives to nominal values
%         dtheta_dn1 = 0;
%         dtheta_dn2 = 0;
% 
%         if n3 > 0  % Near north pole (θ ≈ 0)
%             dtheta_dn3 = -1;  % Approaching from positive side
%         else       % Near south pole (θ ≈ π)
%             dtheta_dn3 = 1;   % Approaching from negative side
%         end
%     else
%         % For components (n1,n2), we have:
%         % ∂n3/∂n_i = (∂/∂n_i)(n3/r) = (∂n3/∂n_i)/r - n3/(r^2) * (∂r/∂n_i)
%         % For n3 component itself:
%         % ∂n3/∂n3 = (∂/∂n3)(n3/r) = 1/r - n3/(r^2) * (n3/r)
% 
%         % Calculate ∂r/∂n_i for each parameter
%         dr_dn1 = n1; % ∂r/∂n1 = n1/r, but r has been normalized
%         dr_dn2 = n2; % ∂r/∂n2 = n2/r, but r has been normalized
%         dr_dn3 = n3; % ∂r/∂n3 = n3/r, but r has been normalized
% 
%         % Now calculate ∂n3/∂n_i
%         dn3_dn1 = -n3 * dr_dn1 / r;
%         dn3_dn2 = -n3 * dr_dn2 / r;
%         dn3_dn3 = 1/r - n3 * dr_dn3 / r;
% 
%         % Finally, calculate ∂θ/∂n_i
%         dtheta_dn1 = -dn3_dn1 / den_theta;
%         dtheta_dn2 = -dn3_dn2 / den_theta;
%         dtheta_dn3 = -dn3_dn3 / den_theta;
%     end
% 
%     % Partial derivatives for phi (atan2(n2, n1))
%     % ∂φ/∂n_i = (∂/∂n_i)atan2(n2, n1)
% 
%     % The derivative of atan2(y, x) is:
%     % ∂/∂x atan2(y, x) = -y/(x²+y²)
%     % ∂/∂y atan2(y, x) = x/(x²+y²)
% 
%     den_phi = n1^2 + n2^2;
% 
%     % Handle numerical issues when n1 and n2 are both near zero
%     if den_phi < 1e-10
%         % When close to the poles (n1≈0, n2≈0), phi becomes ill-defined
%         % Set derivatives to zero since phi changes have minimal effect near poles
%         dphi_dn1 = 0;
%         dphi_dn2 = 0;
%         dphi_dn3 = 0;
%     else
%         % Calculate normalized derivatives of phi
% 
%         % For components n1 and n2 directly involved in atan2:
%         % ∂φ/∂n1 = (∂/∂n1)atan2(n2, n1) = -n2/(n1²+n2²)
%         dphi_dn1_raw = -n2 / den_phi;
% 
%         % ∂φ/∂n2 = (∂/∂n2)atan2(n2, n1) = n1/(n1²+n2²)
%         dphi_dn2_raw = n1 / den_phi;
% 
%         % For n3, changes don't directly affect φ in the atan2 function,
%         % but they affect normalization:
%         dphi_dn3_raw = 0;
% 
%         % Now account for normalization effects on n1, n2
%         dn1_dn1 = 1/r - n1 * dr_dn1 / r;
%         dn1_dn2 = -n1 * dr_dn2 / r;
%         dn1_dn3 = -n1 * dr_dn3 / r;
% 
%         dn2_dn1 = -n2 * dr_dn1 / r;
%         dn2_dn2 = 1/r - n2 * dr_dn2 / r;
%         dn2_dn3 = -n2 * dr_dn3 / r;
% 
%         % Final derivatives using chain rule
%         dphi_dn1 = dphi_dn1_raw * dn1_dn1 + dphi_dn2_raw * dn2_dn1;
%         dphi_dn2 = dphi_dn1_raw * dn1_dn2 + dphi_dn2_raw * dn2_dn2;
%         dphi_dn3 = dphi_dn1_raw * dn1_dn3 + dphi_dn2_raw * dn2_dn3;
%     end
% 
%     % Construct the Jacobian matrix
%     J = [dtheta_dn1, dtheta_dn2, dtheta_dn3;
%          dphi_dn1, dphi_dn2, dphi_dn3];
% end
% 
% 
% function currentFitPSF = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % This function calculates PSF using Cartesian dipole orientation representation
%     % Convert Cartesian to spherical representation for the models that require it
% 
%     % Ensure parameters are within their allowed ranges
%     newangle1 = max(min(newangle1, 1), -1);  % Ensure newangle1 is in [-1,1]
%     newangle2 = max(min(newangle2, 1), -1);  % Ensure newangle2 is in [-1,1]
%     newangle3 = max(min(newangle3, 1), 0);   % Ensure newangle3 is in [0,1]
% 
%     % If it's a Gaussian model, we don't need orientations
%     if strcmpi(model, 'gaussian')
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
%     else
%         % For models that use dipole orientation
%         paramEst.position = Length([x, y, 0], 'nm');
%         paramEst.defocus = Length(defocus, 'nm');
% 
%         % Convert Cartesian coordinates to spherical coordinates (theta, phi)
%         % assuming unit vector normalization
%         % Calculate theta and phi from Cartesian coordinates
%         normFactor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
% 
%         % Normalize to ensure we have a unit vector
%         newangle1_norm = newangle1 / normFactor;
%         newangle2_norm = newangle2 / normFactor;
%         newangle3_norm = newangle3 / normFactor;
% 
%         % Calculate theta (inclination angle from z-axis)
%         theta = acos(newangle3_norm);
% 
%         % Calculate phi (azimuthal angle in xy-plane)
%         phi = atan2(newangle2_norm, newangle1_norm);
% 
%         % Ensure phi is in [0, 2π]
%         if phi < 0
%             phi = phi + 2*pi;
%         end
% 
%         paramEst.dipole = Dipole(theta, phi);
%     end
% 
%     if strcmpi(model, 'mortensen')
%         % Use Mortensen model
%         % Try two specific locations for the Python module (local vs cluster)
%         possibleDirs = {
%             fullfile(pwd, '..', 'simulation-and-fit', 'mortensen_python_bits'), ...
%             '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/mortensen_python_bits' 
%         };
% 
%         % Try each location until we find the Python file
%         pyModuleFound = false;
%         for i = 1:length(possibleDirs)
%             pyDir = possibleDirs{i};
% 
%             % Check if directory exists
%             if ~exist(pyDir, 'dir')
%                 continue;  % Skip to next directory
%             end
% 
%             % Check if Python file exists in this directory
%             pyFilePath = fullfile(pyDir, 'vectorized_mortensen_flipy_phi_theta.py');
%             if exist(pyFilePath, 'file')
%                 % Add to Python path
%                 if count(py.sys.path(), pyDir) == 0
%                     py.sys.path().insert(int32(0), pyDir);
%                 end
% 
%                 pyModuleFound = true;
%                 break;  % Found the file, stop looking
%             end
%         end
% 
%         % If the module wasn't found in any location, show an error
%         if ~pyModuleFound
%             error(['Python module "vectorized_mortensen_flipy_phi_theta.py" not found in any of the specified locations. ', ...
%                    'Please ensure the file exists in one of these directories:\n', ...
%                    '- %s\n', ...
%                    '- %s'], possibleDirs{1}, possibleDirs{2});
%         end
% 
%         % Run the function defined in your python file using the calculated theta and phi
%         currentPsf = py.vectorized_mortensen_flipy_phi_theta.run_simulator_vectorized( ...
%             x, ... % x
%             y, ... % y
%             theta, ... % theta
%             phi, ... % phi
%             paramEst.nPixels,... % image_size_px
%             double(paramEst.pixelSize.inNanometer), ... % pixel_size_nm
%             double(paramEst.wavelength.inNanometer), ... % wavelength
%             paramEst.refractiveIndices(2), ... % n_objective
%             paramEst.refractiveIndices(1), ... % n_sample
%             paramEst.objectiveNA, ... % NA
%             photons ... % n_photons
%         );
% 
%         % Convert to matlab array
%         py_shape = currentPsf.shape;
%         rows = double(py_shape{1});
%         cols = double(py_shape{2});
% 
%         % Initialize MATLAB array
%         psf_matlab = zeros(rows, cols);
% 
%         % Copy values individually using item() method with explicit integer conversion
%         for i = 0:(rows-1)
%             for j = 0:(cols-1)
%                 % Use py.int to explicitly convert indices to Python integers
%                 psf_matlab(i+1, j+1) = double(currentPsf.item(py.tuple({py.int(i), py.int(j)})));
%             end
%         end
% 
%         currentPsf = psf_matlab;
% 
%         % this bit is the equivalent of the getIntensitiesCamera()
%         % stuff done in hinterer
%         totalIntensity = sum(sum(currentPsf));
%         currentPsf = currentPsf / totalIntensity * photons;
%         currentFitPSF = currentPsf;
% 
%     else
%         % Select appropriate BackFocalPlane function based on model
%         if strcmpi(model, 'gaussian')
%             bfp = BackFocalPlane_gaussian(paramEst); % Use Gaussian model
%         elseif strcmpi(model, 'hinterer')
%             bfp = BackFocalPlane(paramEst); % Use Hinterer model
%         end
% 
%         paramEst.backFocalPlane = bfp;
% 
%         % dave
%         paramEst.fieldBFP.x = bfp.electricField.x;
%         paramEst.fieldBFP.y = bfp.electricField.y;
% 
%         currentPsf = zeros(paramEst.nPixels,paramEst.nPixels); 
%         for k=1:size(paramEst.stageDrift.motion,1)
%             % Apply aberrations
%             aberrationCoeffs = getAberrations(paramEst,k);
%             fieldBFP = applyAberrations(paramEst, aberrationCoeffs);
%             % Get image from BFP field
%             currentPsf = currentPsf + getIntensitiesCamera(paramEst, fieldBFP);
%         end
% 
%         totalIntensity = sum(currentPsf,'all');
%         currentPsf = currentPsf ./ totalIntensity * photons + noise_est;
%         currentFitPSF = currentPsf;
%     end
% 
%     currentFitPSF = adjustExcitation(paramEst, currentFitPSF);
%     currentFitPSF = applyShotNoise(paramEst, currentFitPSF);
%     currentFitPSF = addBackgroundNoise(paramEst, currentFitPSF);
% end
% 
% function derivatives = calculateDerivativesCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model)
%     % Calculate derivatives of the PSF with respect to each parameter including newangle1, newangle2, newangle3
%     % Using finite differences for approximation
% 
%     % Small step sizes for numerical derivatives
%     delta_pos = 1;     % nm
%     delta_angle = 0.01; % for Cartesian coordinates
% 
%     % Calculate PSF at current parameters
%     psf0 = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
% 
%     % Derivatives with respect to x and y position
%     psf_dx_plus = calculateProbabilityCartesian(paramEst, x + delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dx_minus = calculateProbabilityCartesian(paramEst, x - delta_pos, y, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dx = (psf_dx_plus - psf_dx_minus) / (2 * delta_pos);
% 
%     psf_dy_plus = calculateProbabilityCartesian(paramEst, x, y + delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     psf_dy_minus = calculateProbabilityCartesian(paramEst, x, y - delta_pos, defocus, newangle1, newangle2, newangle3, photons, noise_est, model);
%     dy = (psf_dy_plus - psf_dy_minus) / (2 * delta_pos);
% 
%     % For models that only use position (like Gaussian)
%     if strcmpi(model, 'gaussian')
%         derivatives = {dx, dy};
%     else
%         % Derivative with respect to newangle1 (with boundary checks)
%         if newangle1 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf0 - psf_da1_minus) / delta_angle;
%         elseif newangle1 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da1_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 + delta_angle, newangle2, newangle3, photons, noise_est, model);
%             psf_da1_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1 - delta_angle, newangle2, newangle3, photons, noise_est, model);
%             da1 = (psf_da1_plus - psf_da1_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle2 (with boundary checks)
%         if newangle2 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf0 - psf_da2_minus) / delta_angle;
%         elseif newangle2 < -1 + delta_angle
%             % Near lower boundary, use forward difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da2_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 + delta_angle, newangle3, photons, noise_est, model);
%             psf_da2_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2 - delta_angle, newangle3, photons, noise_est, model);
%             da2 = (psf_da2_plus - psf_da2_minus) / (2 * delta_angle);
%         end
% 
%         % Derivative with respect to newangle3 (with boundary checks)
%         if newangle3 > 1 - delta_angle
%             % Near upper boundary, use backward difference
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf0 - psf_da3_minus) / delta_angle;
%         elseif newangle3 < delta_angle
%             % Near lower boundary, use forward difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf0) / delta_angle;
%         else
%             % Use central difference
%             psf_da3_plus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 + delta_angle, photons, noise_est, model);
%             psf_da3_minus = calculateProbabilityCartesian(paramEst, x, y, defocus, newangle1, newangle2, newangle3 - delta_angle, photons, noise_est, model);
%             da3 = (psf_da3_plus - psf_da3_minus) / (2 * delta_angle);
%         end
% 
%         % Include orientation derivatives for dipole models
%         derivatives = {dx, dy, da1, da2, da3};
%     end
% end
