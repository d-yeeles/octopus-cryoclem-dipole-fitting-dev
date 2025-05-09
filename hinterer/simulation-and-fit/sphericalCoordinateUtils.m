% function dist = angularDistance(angle1, angle2, period)
%     % Calculate the minimum distance between two angles accounting for periodicity
%     % 
%     % Parameters:
%     % - angle1, angle2: The two angles to compare
%     % - period: The period (2*pi for phi, pi for theta)
%     %
%     % Returns:
%     % - The minimum angular distance
% 
%     % Normalize angles to [0, period)
%     angle1 = mod(angle1, period);
%     angle2 = mod(angle2, period);
% 
%     % Calculate direct distance
%     direct_distance = abs(angle1 - angle2);
% 
%     % Calculate distance going the other way around the circle
%     wrapped_distance = period - direct_distance;
% 
%     % Return the smaller of the two distances
%     dist = min(direct_distance, wrapped_distance);
% end
% 
% function [theta, phi] = cartesianToSpherical(newangle1, newangle2, newangle3)
%     % Convert from Cartesian unit vector to spherical coordinates
% 
%     % Ensure the vector is normalized
%     norm_factor = sqrt(newangle1^2 + newangle2^2 + newangle3^2);
%     if norm_factor > 1e-10
%         newangle1 = newangle1 / norm_factor;
%         newangle2 = newangle2 / norm_factor;
%         newangle3 = newangle3 / norm_factor;
%     end
% 
%     % Calculate spherical coordinates
%     theta = acos(newangle3);                  % theta in [0, π]
%     phi = atan2(newangle2, newangle1);        % phi in [-π, π]
%     phi = mod(phi, 2*pi);                     % phi in [0, 2π]
% end
% 
% function [newangle1, newangle2, newangle3] = sphericalToCartesian(theta, phi)
%     % Convert from spherical coordinates to Cartesian unit vector
% 
%     % Ensure angles are in the correct ranges
%     theta = mod(theta, pi);                  % theta in [0, π]
%     phi = mod(phi, 2*pi);                    % phi in [0, 2π]
% 
%     % Calculate Cartesian coordinates
%     newangle1 = sin(theta) * cos(phi);       % x-component
%     newangle2 = sin(theta) * sin(phi);       % y-component
%     newangle3 = cos(theta);                  % z-component
% end
% 
% function angular_covariance = transformCovarianceToAngular(covariance, newangle1, newangle2, newangle3)
%     % Transform the covariance matrix from Cartesian to spherical coordinates
% 
%     % Extract position and orientation parts
%     position_cov = covariance(1:2, 1:2);
%     orientation_cov = covariance(3:5, 3:5);
%     position_orientation_cov = covariance(1:2, 3:5);
% 
%     % Normalize the orientation vector
%     n = [newangle1, newangle2, newangle3];
%     norm_value = norm(n);
%     if norm_value < 1e-10
%         error('Orientation vector is too close to zero');
%     end
%     n = n / norm_value;
%     newangle1 = n(1);
%     newangle2 = n(2);
%     newangle3 = n(3);
% 
%     % Calculate Jacobian for the transformation
%     % This is a 2x3 matrix of partial derivatives of [theta, phi] with respect to [x, y, z]
% 
%     % Precompute some values
%     xy_squared = newangle1^2 + newangle2^2;
%     r_squared = xy_squared + newangle3^2;
%     xy_norm = sqrt(xy_squared);
% 
%     % Avoid division by zero near the poles
%     if xy_norm < 1e-10
%         xy_norm = 1e-10;
%     end
% 
%     % Partial derivatives for theta
%     dtheta_dx = -newangle1 * newangle3 / (r_squared * xy_norm);
%     dtheta_dy = -newangle2 * newangle3 / (r_squared * xy_norm);
%     dtheta_dz = xy_norm / r_squared;
% 
%     % Partial derivatives for phi
%     dphi_dx = -newangle2 / xy_squared;
%     dphi_dy = newangle1 / xy_squared;
%     dphi_dz = 0;  % phi doesn't depend on z
% 
%     % Construct the Jacobian matrix
%     J = [dtheta_dx, dtheta_dy, dtheta_dz;
%          dphi_dx, dphi_dy, dphi_dz];
% 
%     % Transform orientation covariance using the Jacobian
%     angular_orientation_cov = J * orientation_cov * J';
% 
%     % Transform the position-orientation cross-covariance
%     position_angular_cov = position_orientation_cov * J';
% 
%     % Construct the full transformed covariance matrix
%     angular_covariance = [position_cov, position_angular_cov;
%                          position_angular_cov', angular_orientation_cov];
% end