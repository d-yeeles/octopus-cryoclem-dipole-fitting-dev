classdef BackFocalPlane_gaussian

    % Calculates the BFP-fields (Ex,Ey) of an arbitrary oriented dipole in
    % the focal point of an objective based on the following paper:
    % Axelrod, J. of Microscopy 2012
    %
    % Details:
    % - An intermediate layer between coverglass and specimen can be assumed
    % - Includes also near-field effects (SAF emission)

    properties
        nGrid (1,1) {mustBeInteger, mustBePositive} = 129
        unitKSpace (1,1) {mustBePositive} = 1
        electricField
    end

    methods
        function obj = BackFocalPlane_gaussian(psfObj)
            obj.nGrid = psfObj.nDiscretizationBFP;
            obj.unitKSpace = psfObj.unitKSpace; % k-space unit in BFP
            obj.electricField = calculateElectricField(obj,psfObj);
        end

        function E_BFP = calculateElectricField(obj, psfObj)
        

            % This version probably wrong

            % Coordinates in the objective pupil
            par.nGrid = obj.nGrid;
            par.spacing = obj.unitKSpace;
            pupilMask = Mask(par);
            [Kr, ~] = pupilMask.getPolarCoordinates; % Normalised radial coordinates in the pupil
            % Kr = pupilMask.radius; % Not normalised

            % Define Gaussian parameters
            sigma = 0.3; % Standard deviation of the Gaussian
            amplitude = 1; % Peak amplitude of the Gaussian

            % Gaussian Electric Field
            % GaussianField = abs(amplitude * exp(-Kr.^2 / (2 * sigma^2)));
            GaussianField = amplitude * exp(-Kr.^2 / (2 * sigma^2));

            % Output Gaussian electric field for both x and y components
            E_BFP.x = GaussianField; % x-component of the electric field
            E_BFP.y = GaussianField; % y-component of the electric field (same as x for symmetry)
        





            % % This might be more correct (amplitude and sigma not just set
            % % to 1 and 0.5):
            % % Work it out fully, then take the max amplitude of that
            % 
            % RI = psfObj.refractiveIndices;
            % dipole = psfObj.dipole;
            % hIntermediate = psfObj.heightIntermediateLayer.inMeter;
            % pos = psfObj.position.inMeter;
            % z = pos(3);
            % focalLength = psfObj.objectiveFocalLength.inMeter;
            % mu = 1e-12;
            % 
            % %% Pre-Calculations
            % 
            % if length(RI)==1
            %     RI=[RI, RI, RI];
            % end
            % 
            % % Coordinates in the objective pupil
            % par.nGrid = obj.nGrid;
            % par.spacing = obj.unitKSpace;
            % pupilMask = Mask(par);
            % [~, maskAngle] = pupilMask.getPolarCoordinates;
            % PHI3 = fliplr(maskAngle) - pi;
            % 
            % % Wavenumbers (magnitude of k-vectors) in the different media, k = k0 * RI
            % k0 = 2 * pi / psfObj.wavelength.inMeter; % wavenumber in vacuum
            % k1 = k0 * RI(1); % wavenumber in media 1 (typically water), Eq. (14) of paper
            % k2 = k0 * RI(2); % wavenumber in media 2 (intermediate layer), Eq. (14) of paper
            % k3 = k0 * RI(3); % wavenumber in media 3 (immersion medium), Eq. (14) of paper
            % 
            % % Angles in different media
            % Kr = pupilMask.radius;
            % THETA1 = acos( sqrt( 1 - (RI(3)/RI(1) * Kr/k3).^2 ) ) .* pupilMask.values; % angle in medium 1; Eq. (4) of paper
            % THETA2 = acos( sqrt( 1 - (RI(3)/RI(2) * Kr/k3).^2 ) ) .* pupilMask.values; % angle in medium 2; Eq. (4) of paper
            % THETA3 = asin( Kr/k3 ) .* pupilMask.values; % angle in medium 3; maximum theta3 from Eq. (19) of paper
            % 
            % %% Calculations according to paper of Axelrod, 2012
            % 
            % % Cosines of angles
            % CTHETA1 = cos(THETA1);
            % CTHETA2 = cos(THETA2);
            % CTHETA3 = cos(THETA3);
            % 
            % % Fresnel-coefficients
            % % Eq. (3) of paper
            % tp12 = 2*RI(1)*CTHETA1./(RI(1)*CTHETA2+RI(2)*CTHETA1);
            % tp23 = 2*RI(2)*CTHETA2./(RI(2)*CTHETA3+RI(3)*CTHETA2);
            % 
            % ts12 = 2*RI(1)*CTHETA1./(RI(1)*CTHETA1+RI(2)*CTHETA2);
            % ts23 = 2*RI(2)*CTHETA2./(RI(2)*CTHETA2+RI(3)*CTHETA3);
            % 
            % rp12 = (RI(2)*CTHETA1-RI(1)*CTHETA2)./(RI(1)*CTHETA2+RI(2)*CTHETA1);
            % rp23 = (RI(3)*CTHETA2-RI(2)*CTHETA3)./(RI(2)*CTHETA3+RI(3)*CTHETA2);
            % 
            % rs12 = (RI(1)*CTHETA1-RI(2)*CTHETA2)./(RI(1)*CTHETA1+RI(2)*CTHETA2);
            % rs23 = (RI(2)*CTHETA2-RI(3)*CTHETA3)./(RI(2)*CTHETA2+RI(3)*CTHETA3);
            % 
            % % Fresnel coefficients for three-layer system
            % % Eq. (12) of paper
            % tp = tp12 .* tp23 .* exp(1i*k2*hIntermediate*CTHETA2) ./ (1 + rp12 .* rp23 .* exp(2i*k2*hIntermediate*CTHETA2));
            % ts = ts12 .* ts23 .* exp(1i*k2*hIntermediate*CTHETA2) ./ (1 + rs12 .* rs23 .* exp(2i*k2*hIntermediate*CTHETA2));
            % 
            % % Dipole projections onto directions p, s and z
            % % Eq. (13) of paper
            % mu_p = mu * sin(dipole.inclination) .* cos(dipole.azimuth - PHI3);
            % mu_s = mu * sin(dipole.inclination) .* sin(dipole.azimuth - PHI3);
            % mu_z = mu * cos(dipole.inclination);
            % 
            % % Prefactor C (the only constant where f plays a role)
            % % Eq. (11) of paper
            % C = ( k3^2 * exp(1i*k3*focalLength) .* CTHETA3) / (focalLength * RI(1)) ...
            %     .* exp(-1i*k3*hIntermediate*CTHETA3) ...
            %     .* exp(1i*k1.*CTHETA1.*z);
            % 
            % % Electric field components in layer 3 (pre-objective zone), along the p, s and z-axes
            % % (see paper for axes definitions)
            % % Eq. (10) of paper
            % E3p = C .* tp .* CTHETA3 .* (mu_p./RI(3) + mu_z.*sin(THETA3)./CTHETA1);
            % E3s = C .* ts .* (mu_s./(RI(3)./CTHETA1));
            % E3z = C .* tp .* sin(THETA3) .* (mu_p./RI(3) + mu_z.*sin(THETA3)./CTHETA1);
            % 
            % % Influence of objective, rotation of rays by their angle theta3 such that they are all parallel to the optical axis
            % % Eq. (15) and (16) of paper
            % apodizationFactor = 1 ./ sqrt(CTHETA3) .* pupilMask.values; % apodization of objective lens
            % E_BFP_p = (E3p.*CTHETA3 + E3z.*sin(THETA3)) .* apodizationFactor;
            % E_BFP_s = E3s .* apodizationFactor;  % s-polarization remains unchanged by this rotation
            % 
            % % Coordinate transformation into x-and y-polarization yields fields in the back focal plane of objective
            % % Eq. (18) of paper
            % E_BFP.x = cos(PHI3).*E_BFP_p - sin(PHI3).*E_BFP_s;
            % E_BFP.y = sin(PHI3).*E_BFP_p + cos(PHI3).*E_BFP_s;
            % 
            % % Now construct Gaussian using the above to derive peak
            % % amplitude and std
            % I_total = abs(E_BFP.x).^2 + abs(E_BFP.y).^2;
            % I_max = max(max(I_total));
            % half_I_max = I_max / 2;
            % % Identify the radial distances where the intensity is less than or equal to half of the maximum
            % intensity_below_half = I_total <= half_I_max;
            % % Extract the radial distance corresponding to this intensity (assuming Kr corresponds to radial distance)
            % Kr_at_FWHM = Kr(intensity_below_half);  % Kr is the radial coordinate
            % % Calculate the FWHM (distance between points where intensity is half maximum)
            % FWHM = max(Kr_at_FWHM) - min(Kr_at_FWHM);
            % % Calculate sigma from FWHM
            % sigma = FWHM / 2.355;
            % amplitude = I_max;
            % 
            % % Gaussian Electric Field - restrict to positive values
            % GaussianField = amplitude * exp(-Kr.^2 / (2 * sigma^2));
            % 
            % % Output Gaussian electric field for both x and y components
            % E_BFP.x = GaussianField; % x-component of the electric field
            % E_BFP.y = GaussianField; % y-component of the electric field (same as x for symmetry)
            % 
            % % disp(amplitude)
            % % disp(sigma)


            % % Extract parameters from the psfObj
            % RI = psfObj.refractiveIndices;
            % pos = psfObj.position.inMeter; % position in meters
            % focalLength = psfObj.objectiveFocalLength.inMeter;
            % wavelength = psfObj.wavelength.inMeter;
            % 
            % %% Ensure RI has three values
            % if length(RI) == 1
            %     RI = [RI, RI, RI];
            % end
            % 
            % %% Pre-Calculations
            % % Coordinates in the objective pupil
            % par.nGrid = obj.nGrid;
            % par.spacing = obj.unitKSpace;
            % pupilMask = Mask(par);
            % Kr = pupilMask.radius;
            % 
            % % Define Gaussian parameters
            % gaussianCenterPolar = pos(1:2); % center in polar coordinates (r, theta)
            % gaussianSigma = 400 / 2.355; % standard deviation (adjustable parameter). I think this needs to be in metres because pos is in metres
            % 
            % % Shift polar coordinates for Gaussian center
            % deltaKr = Kr - gaussianCenterPolar(1); % radial distance shift
            % 
            % % Compute the Gaussian field in the back focal plane
            % gaussianField = exp(-(deltaKr.^2) / (2 * gaussianSigma^2));
            % 
            % % Combine Gaussian envelope with phase factor
            % E_BFP.x = gaussianField;
            % E_BFP.y = gaussianField;




            % % Extract parameters from the psfObj
            % pos = psfObj.position.inMeter; % position in meters
            % focalLength = psfObj.objectiveFocalLength.inMeter;
            % 
            % %% Pre-Calculations
            % % Coordinates in the objective pupil
            % par.nGrid = obj.nGrid;
            % par.spacing = obj.unitKSpace;
            % pupilMask = Mask(par);
            % [~, maskAngle] = pupilMask.getPolarCoordinates();
            % Kr = pupilMask.radius;
            % 
            % % Define Gaussian parameters
            % gaussianCenter = pos(1:2); % center in x and y from pos
            % gaussianSigma = (400e30)/(2*sqrt(2*log(2))); % standard deviation (FWHM of a blob, guessed at 400nm. think this should be in m though.)
            % 
            % % Convert polar coordinates to Cartesian for Gaussian calculation
            % [X, Y] = pol2cart(maskAngle, Kr);
            % 
            % % Compute the Gaussian field in the back focal plane
            % gaussianEnvelope = exp(-((X - gaussianCenter(1)).^2 + (Y - gaussianCenter(2)).^2) / (2 * gaussianSigma^2));
            % 
            % % Combine Gaussian envelope with phase factor
            % E_BFP.x = gaussianEnvelope;
            % E_BFP.y = gaussianEnvelope;




        end


        function plot(obj)
            %% Electric field
            % Electric field x
            subplot(2,2,1)
            imagesc(real(obj.electricField.x))
            axis equal; axis tight
            title('E_x real part')

            subplot(2,2,2)
            imagesc(imag(obj.electricField.x))
            axis equal; axis tight
            title('E_x imaginary part')

            % Electric field y
            subplot(2,2,3)
            imagesc(real(obj.electricField.y))
            axis equal; axis tight
            title('E_y real part')

            subplot(2,2,4)
            imagesc(imag(obj.electricField.y))
            axis equal; axis tight
            title('E_y imaginary part')

            %% Intensity
            figure
            subplot(1,3,1)
            imagesc(abs(obj.electricField.x).^2)
            axis equal; axis tight
            title('Intensity x-pol.')

            subplot(1,3,2)
            imagesc(abs(obj.electricField.y).^2)
            axis equal; axis tight
            title('Intensity y-pol.')

            subplot(1,3,3)
            imagesc(abs(obj.electricField.x).^2 + abs(obj.electricField.y).^2)
            axis equal; axis tight
            title('Total intensity')
        end
    end
end