classdef Aperture < Attenuation
    
    methods
        function obj = Aperture(gridSize)
            sector = ConstantAttenuationMask(gridSize,0);
            obj.mask = sector.mask;
        end
    end
endclassdef Attenuation

    properties
        mask (:,:) double {mustBeInRange(mask,0,1)}
    end

    properties (Dependent, Hidden)
        gridSize
        polarCoordinates
    end


    methods
        %% Constructor
        function obj = Attenuation(mask)
            if nargin > 0
                obj.mask = mask;
            end
        end

        %% Get methods
        function mask = getMask(obj)
            mask = obj.mask;
        end
        function gridSize = get.gridSize(obj)
            gridSize = size(obj.mask,1);
        end
        function polarCoord = get.polarCoordinates(obj)
            [theta,rho] = Attenuation.createGrid(obj.gridSize);
            polarCoord.angle = theta;
            polarCoord.radius = rho;
        end
    end

    methods (Access=public)
        %% Functions
        function attenuatedElectricField = apply(obj, electricField)
            assert(isequal(size(obj.mask),size(electricField)))
            attenuatedElectricField = sqrt(1-obj.mask) .* electricField;
            % Note: Full attenuation (=1) ~= no electric field (.*0)
            % Attenuation coefficient is for intensity;
            % take sqrt to apply attenuation directly to electric field
        end

        %% Plot
        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            imagesc(axesHandle, obj.mask)
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            colormap(axesHandle, flipud(gray))
            caxis(axesHandle, [0 1])
            cb.Label.String = 'Attenuation';
            cb.Label.FontSize = 12;
            cb.Ticks = (0:0.25:1);
            cb.TickLabels = {'0','','0.5','','1'};
            set(axesHandle,'visible','off')
        end
    end

    methods (Access=public)
        %% Adapt attenuation mask
        function obj = cutInnerRing(obj, relativeRadius)
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius < relativeRadius) = 1;
        end

        function obj = cutAperture(obj, relativeRadius)
            if nargin < 2
                relativeRadius = 1;
            end
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius > relativeRadius) = 1;
        end

        function obj = selectCircleSector(obj, sectorAngles, attenuation)
            mustBeInFullRadialRange(sectorAngles)
            if numel(sectorAngles)==1
                sectorAngles = [0,sectorAngles];
            end
            if nargin < 3
                attenuation = 1;
            end
            isInSegment = ((sectorAngles(1)<obj.polarCoordinates.angle)) & (obj.polarCoordinates.angle<sectorAngles(2));
            obj.mask(isInSegment) = attenuation;
            obj = obj.cutAperture();
        end

        function obj = rotate(obj,rotationAngle)
            mustBeInFullRadialRange(rotationAngle)
            if rotationAngle ~= 0
                obj.mask = imrotate(obj.mask, -rotationAngle*180/pi, 'bilinear', 'crop');
                % Fill cropped area with ones
                obj.mask(obj.polarCoordinates.radius > 1) = 1;
            end
        end

        %% Overload operators
        function combinedPhaseMask = plus(obj1,obj2)
            combinedPhaseMask = PhaseMask(obj1.mask + obj2.mask);
        end

        function combinedPhaseMask = times(obj1,obj2)
            combinedPhaseMask = PhaseMask(obj1.mask .* obj2.mask);
        end
    end


    methods (Static)
        function [theta,rho] = createGrid(gridSize)
            x = linspace(-1, 1, gridSize);
            [theta, rho] = cart2pol(x', x);
            theta = flipud(mod(theta,2*pi));
        end
    end
end
classdef BackFocalPlane

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
        function obj = BackFocalPlane(psfObj)
            obj.nGrid = psfObj.nDiscretizationBFP;
            obj.unitKSpace = psfObj.unitKSpace; % k-space unit in BFP
            obj.electricField = calculateElectricField(obj,psfObj);
        end

        function E_BFP = calculateElectricField(obj,psfObj)

            RI = psfObj.refractiveIndices;
            dipole = psfObj.dipole;
            hIntermediate = psfObj.heightIntermediateLayer.inMeter;
            pos = psfObj.position.inMeter;
            z = pos(3);
            focalLength = psfObj.objectiveFocalLength.inMeter;
            mu = 1e-12;

            %% Pre-Calculations

            if length(RI)==1
                RI=[RI, RI, RI];
            end

            % Coordinates in the objective pupil
            par.nGrid = obj.nGrid;
            par.spacing = obj.unitKSpace;
            pupilMask = Mask(par);
            [~, maskAngle] = pupilMask.getPolarCoordinates;
            PHI3 = fliplr(maskAngle) - pi;

            % Wavenumbers (magnitude of k-vectors) in the different media, k = k0 * RI
            k0 = 2 * pi / psfObj.wavelength.inMeter; % wavenumber in vacuum
            k1 = k0 * RI(1); % wavenumber in media 1 (typically water), Eq. (14) of paper
            k2 = k0 * RI(2); % wavenumber in media 2 (intermediate layer), Eq. (14) of paper
            k3 = k0 * RI(3); % wavenumber in media 3 (immersion medium), Eq. (14) of paper

            % Angles in different media
            Kr = pupilMask.radius;
            THETA1 = acos( sqrt( 1 - (RI(3)/RI(1) * Kr/k3).^2 ) ) .* pupilMask.values; % angle in medium 1; Eq. (4) of paper
            THETA2 = acos( sqrt( 1 - (RI(3)/RI(2) * Kr/k3).^2 ) ) .* pupilMask.values; % angle in medium 2; Eq. (4) of paper
            THETA3 = asin( Kr/k3 ) .* pupilMask.values; % angle in medium 3; maximum theta3 from Eq. (19) of paper

            %% Calculations according to paper of Axelrod, 2012

            % Cosines of angles
            CTHETA1 = cos(THETA1);
            CTHETA2 = cos(THETA2);
            CTHETA3 = cos(THETA3);

            % Fresnel-coefficients
            % Eq. (3) of paper
            tp12 = 2*RI(1)*CTHETA1./(RI(1)*CTHETA2+RI(2)*CTHETA1);
            tp23 = 2*RI(2)*CTHETA2./(RI(2)*CTHETA3+RI(3)*CTHETA2);

            ts12 = 2*RI(1)*CTHETA1./(RI(1)*CTHETA1+RI(2)*CTHETA2);
            ts23 = 2*RI(2)*CTHETA2./(RI(2)*CTHETA2+RI(3)*CTHETA3);

            rp12 = (RI(2)*CTHETA1-RI(1)*CTHETA2)./(RI(1)*CTHETA2+RI(2)*CTHETA1);
            rp23 = (RI(3)*CTHETA2-RI(2)*CTHETA3)./(RI(2)*CTHETA3+RI(3)*CTHETA2);

            rs12 = (RI(1)*CTHETA1-RI(2)*CTHETA2)./(RI(1)*CTHETA1+RI(2)*CTHETA2);
            rs23 = (RI(2)*CTHETA2-RI(3)*CTHETA3)./(RI(2)*CTHETA2+RI(3)*CTHETA3);

            % Fresnel coefficients for three-layer system
            % Eq. (12) of paper
            tp = tp12 .* tp23 .* exp(1i*k2*hIntermediate*CTHETA2) ./ (1 + rp12 .* rp23 .* exp(2i*k2*hIntermediate*CTHETA2));
            ts = ts12 .* ts23 .* exp(1i*k2*hIntermediate*CTHETA2) ./ (1 + rs12 .* rs23 .* exp(2i*k2*hIntermediate*CTHETA2));

            % Dipole projections onto directions p, s and z
            % Eq. (13) of paper
            mu_p = mu * sin(dipole.inclination) .* cos(dipole.azimuth - PHI3);
            mu_s = mu * sin(dipole.inclination) .* sin(dipole.azimuth - PHI3);
            mu_z = mu * cos(dipole.inclination);

            % Prefactor C (the only constant where f plays a role)
            % Eq. (11) of paper
            C = ( k3^2 * exp(1i*k3*focalLength) .* CTHETA3) / (focalLength * RI(1)) ...
                .* exp(-1i*k3*hIntermediate*CTHETA3) ...
                .* exp(1i*k1.*CTHETA1.*z);

            % Electric field components in layer 3 (pre-objective zone), along the p, s and z-axes
            % (see paper for axes definitions)
            % Eq. (10) of paper
            E3p = C .* tp .* CTHETA3 .* (mu_p./RI(3) + mu_z.*sin(THETA3)./CTHETA1);
            E3s = C .* ts .* (mu_s./(RI(3)./CTHETA1));
            E3z = C .* tp .* sin(THETA3) .* (mu_p./RI(3) + mu_z.*sin(THETA3)./CTHETA1);

            % Influence of objective, rotation of rays by their angle theta3 such that they are all parallel to the optical axis
            % Eq. (15) and (16) of paper
            apodizationFactor = 1 ./ sqrt(CTHETA3) .* pupilMask.values; % apodization of objective lens
            E_BFP_p = (E3p.*CTHETA3 + E3z.*sin(THETA3)) .* apodizationFactor;
            E_BFP_s = E3s .* apodizationFactor;  % s-polarization remains unchanged by this rotation

            % Coordinate transformation into x-and y-polarization yields fields in the back focal plane of objective
            % Eq. (18) of paper
            E_BFP.x = cos(PHI3).*E_BFP_p - sin(PHI3).*E_BFP_s;
            E_BFP.y = sin(PHI3).*E_BFP_p + cos(PHI3).*E_BFP_s;

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
endclassdef ChirpZTransform
    % 2D chirp-z transform
    %
    % See the following publication:
    % Raoqiong Tong and Robert W. Cox.
    % "Rotation of NMR images using the 2D chirp-z transform"
    % Magnetic resonance in medicine 41.2 (1999): 253-256.
    
    properties
        kernel
        fourierKernel
        nPixelsPadded
    end

    methods
        function obj = ChirpZTransform(psfObj)
            % Get dimensions of image and k-space
            nPixelsImage = psfObj.nPixels .* psfObj.oversampling; % size of field in pixel
            nPixelsBFP = psfObj.nDiscretizationBFP; % size of input field (k-space)
            obj.nPixelsPadded = nPixelsImage + nPixelsBFP - 1; % padded size required for convolution is Nx+Nk-1

            % Create grid
            x = ( floor(-obj.nPixelsPadded(1)/2 + 0.5) : floor(obj.nPixelsPadded(1)/2 - 0.5) )';
            y = ( floor(-obj.nPixelsPadded(1)/2 + 0.5) : floor(obj.nPixelsPadded(1)/2 - 0.5) );

            % Create kernel (quadratic convolution phase kernel)
            alpha = psfObj.unitObjectSpace * psfObj.unitKSpace / (2*pi);
            obj.kernel = exp(-1i*alpha*pi*x.^2) * exp(-1i*alpha*pi*y.^2); % = exp(-1i*alpha*pi*(x.^2+y.^2)), but slightly faster
            
            % Fourier transform of kernel
            obj.fourierKernel = fft2( ifftshift(obj.kernel) );
        end

        function E_out = apply(obj,psfObj,E_in)
            % Complex conjugate of kernel
            conjKernel = conj(obj.kernel);
            
            % Electric field times complex conjugate of kernel
            f = embedArray2D(E_in,obj.nPixelsPadded,0) .* conjKernel;

            % Fourier transform of kernel
            F = fft2( ifftshift(f) );

            % Convolving f with kernel (inverese fft of multiplication)
            convolution = fftshift( ifft2( F .* obj.fourierKernel ) );

            % Final multiplication by factor
            E_out = obj.nPixelsPadded * conjKernel .* convolution;

            % Crop to image size
            E_out = cropArray2D(E_out, psfObj.nPixels .* psfObj.oversampling);
        end
    end
endclassdef ConstantAttenuationMask < Attenuation
    
    methods
        function obj = ConstantAttenuationMask(gridSize, attenuation)
            obj.mask = attenuation*ones(gridSize);
            obj = obj.cutAperture();
        end
    end
endclassdef Dipole

    properties
        inclination (1,1) {mustBeInFullRadialRange} = 0 % Inclination angle
        azimuth (1,1) {mustBeInFullRadialRange} = 0 % Azimuthal angle
    end
    
    methods
        function obj = Dipole(angleInclination, angleAzimuth)
            obj.inclination = angleInclination;
            obj.azimuth = angleAzimuth;
        end
        
        function dipoleVector = getDipoleVector(obj)
            dipoleVector = [obj.inclination, obj.azimuth];
        end
    end
endclassdef EmptyPhaseMask < PhaseMask

    methods
        function obj = EmptyPhaseMask(gridSize)
            obj.mask = zeros(gridSize);
        end
    end

    methods (Access=public)
        % Overwrite apply method:
        function shiftedElectricField = apply(obj, electricField)
            shiftedElectricField = electricField;
        end
    end
end

function [psf, Ix, Iy] = getIntensitiesCamera(obj, fieldBFP)

    Ix = abs( obj.chirpZTransform.apply(obj, fieldBFP.x) ).^2; % intensity = |E_imagePlane|²
    Iy = abs( obj.chirpZTransform.apply(obj, fieldBFP.y) ).^2; % intensity = |E_imagePlane|²

    psf = Ix + Iy;

    % Account for pixel sensitivity
    psf = psf .* repmat(obj.pixelSensitivityMask, obj.nPixels);

    % Sum over block matrices (to return to desired pixelsize)
    psf = squeeze(sum(sum(reshape(psf,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels),1),3));
    
    % Normalization
    totalIntensity = sum(sum(psf));
    psf = psf / totalIntensity * obj.nPhotons;
    Ix = Ix ./ totalIntensity * obj.nPhotons;
    Iy = Iy ./ totalIntensity * obj.nPhotons;
endclassdef Mask

    properties
        values double
        nGrid = 129
        spacing = 1 % spacing between two neighbouring entries
        mode = 'FFT' % 'exact' or 'FFT'
                     % 'exact':
                     % zero point is exactly in the middle,
                     % i.e. for even grid sizes between the two middle pixels
                     % 'FFT':
                     % zero point is the middle pixel for odd grid sizes
                     % and the pixel to the lower right of the exact centre for even grid sizes
                     % (this is also the pixel which MATLAB takes as zero)
        radius
    end

    methods
        function obj = Mask(par)
            if nargin > 0
                obj = setInputParameters('Mask', obj, par);
            end

            % Calculate mask values
            N = [obj.nGrid, obj.nGrid];
            s = [obj.spacing, obj.spacing];
            if strcmp(obj.mode,'exact')==1
                x = linspace( -(N(1)-1)/2, (N(1)-1)/2, N(1) ) * s(1);
                y = linspace( -(N(2)-1)/2, (N(2)-1)/2, N(2) ) * s(2);
            elseif strcmp(obj.mode,'FFT')==1
                x = ( floor( -N(1)/2 + 0.5) : floor( N(1)/2 - 0.5) ) * s(1);
                y = ( floor( -N(2)/2 + 0.5) : floor( N(2)/2 - 0.5) ) * s(2);
            end

            [X, Y] = meshgrid(y,x);
            obj.radius = sqrt(X.^2+Y.^2);
            obj.values = ( obj.radius.^2 <= ((min(N/2.*s))).^2 );
        end

        function [normalizedRadius, angle] = getPolarCoordinates(obj)
            x = linspace(-1,1,obj.nGrid);
            [X,Y] = meshgrid(x,x);
            [angle,normalizedRadius] = cart2pol(X,Y);
            angle = fliplr(angle + pi);

            %idxNan = (obj.values==0);
            %angle(idxNan) = 0; %NaN;
            %radius(idxNan) = 0; %NaN;
        end

        function plot(obj)
            imagesc(obj.values)
            axis equal; axis tight;
        end
    end
endclassdef NoStageDrift < StageDrift

    methods
        function obj = NoStageDrift()
            obj.motion = Length([0 0 0],'nm');
        end
    end
endclassdef PhaseMask

    properties
        mask (:,:) double
    end

    properties (Dependent, Hidden)
        gridSize
        polarCoordinates
    end


    methods
        %% Constructor
        function obj = PhaseMask(mask)
            if nargin > 0
                obj.mask = mask;
            end
        end

        %% Get methods
        function mask = getMask(obj)
            mask = obj.mask;
        end
        function gridSize = get.gridSize(obj)
            gridSize = size(obj.mask,1);
        end
        function polarCoord = get.polarCoordinates(obj)
            [theta,rho] = PhaseMask.createGrid(obj.gridSize);
            polarCoord.angle = theta;
            polarCoord.radius = rho;
        end
    end

    methods (Access=public)
        %% Functions
        function shiftedElectricField = apply(obj, electricField)
            assert(isequal(size(obj.mask),size(electricField)))
            shiftedElectricField = electricField .* exp(1i.*obj.mask);
        end

        %% Plot
        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            imagesc(axesHandle, obj.mask)
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            caxis(axesHandle, [0 2*pi])
            cb.Label.String = 'Phase shift';
            cb.Label.FontSize = 12;
            cb.Ticks = (0:0.5:2)*pi;
            cb.TickLabels = {'0','','\pi','','2\pi'};
            set(axesHandle,'visible','off')
        end
    end

    methods (Access=public)
        %% Adapt phase mask
        function obj = cutInnerRing(obj, relativeRadius)
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius < relativeRadius) = 0;
        end

        function obj = cutAperture(obj, relativeRadius)
            if nargin < 2
                relativeRadius = 1;
            end
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius > relativeRadius) = 0;
        end

        function obj = selectCircleSector(obj, sectorAngles)
            mustBeInFullRadialRange(sectorAngles)
            if numel(sectorAngles)==1
                sectorAngles = [0,sectorAngles];
            end
            isInSegment = ((sectorAngles(1)<obj.polarCoordinates.angle)) & (obj.polarCoordinates.angle<sectorAngles(2));
            obj.mask(~isInSegment) = 0;
        end

        function obj = rotate(obj,rotationAngle)
            mustBeInFullRadialRange(rotationAngle)
            if rotationAngle ~= 0
                obj.mask = imrotate(obj.mask, -rotationAngle*180/pi, 'bilinear', 'crop');
            end
        end

        %% Overload operators
        function combinedPhaseMask = plus(obj1,obj2)
            combinedPhaseMask = PhaseMask(obj1.mask + obj2.mask);
        end

        function combinedPhaseMask = times(obj1,obj2)
            combinedPhaseMask = PhaseMask(obj1.mask .* obj2.mask);
        end
    end

    methods (Static)
        function [theta,rho] = createGrid(gridSize)
            x = linspace(-1, 1, gridSize);
            [theta, rho] = cart2pol(x', x);
            theta = flipud(mod(theta,2*pi));
        end
    end
end
classdef PixelSensitivity
    
    properties
        mask
    end
    
    enumeration
        centerLarge (PixelSensitivity.centerMask(9,1))
        centerSmall (PixelSensitivity.centerMask(9,3))
        centerOnly (PixelSensitivity.centerMask(3,1))
        none (0)
    end
    
    methods
        function obj = PixelSensitivity(mask)
            obj.mask = mask;
        end
    end
    
    methods (Static)
        function mask = uniform(n)
            mask = ones(n);
        end
        
        function mask = centerMask(width, border)
            mask = ones(width);
            mask(:,1:border) = 0;
            mask(:,end-border+1:end) = 0;
            mask(1:border,:) = 0;
            mask(end-border+1:end,:) = 0;
        end
    end
end

classdef PSF

    properties
        % PSF image
        image
        nPixels (1,1) {mustBeInteger, mustBePositive, mustBeOdd} = 17

        % Fluorophore
        dipole Dipole = Dipole(0,0)
        position (1,3) Length = Length([0 0 0], 'nm')
        nPhotons (1,1) {mustBeInteger, mustBeNonnegative} = 1e5
        shotNoise (1,1) logical = 0
        reducedExcitation (1,1) logical = 0
        stageDrift StageDrift = NoStageDrift() % relative lateral motion from start point (3rd component = defocus) 
        
        % Microscope setup
        wavelength (1,1) Length {mustBePositive} = Length(680,'nm')
        defocus (1,1) Length = Length(0, 'nm') % defocus of objective lens (negative values = moving focus into fluid)
        astigmatism = 0 % Zernike coefficient (in units of wavelength, i.e. 0.11 for Zernike coefficient 0.11*lambda)
        objectiveNA (1,1) {mustBePositive} = 0.7
        objectiveFocalLength (1,1) Length {mustBePositive} = Length(3,'mm')
        refractiveIndices (1,3) {mustBePositive} = [1.33 1.46 1] % [RI_specimen, RI_intermed, RI_immoil]
        heightIntermediateLayer (1,1) Length {mustBeNonnegative} = Length(0, 'mm')

        % Back focal plane
        backFocalPlane
        phaseMask = @(n) EmptyPhaseMask(n)
        attenuation = @(n) Aperture(n)
        nDiscretizationBFP (1,1) {mustBeInteger, mustBePositive, mustBeOdd} = 129
        
        % Camera
        pixelSize (1,1) Length {mustBePositive} = Length(108,'nm')
        pixelSensitivityMask {mustBeNonnegative, mustBeOddSize} = PixelSensitivity.uniform(9)
        backgroundNoise (1,1) {mustBeNonnegative} = 0
    end

    properties (Hidden)
        phaseMaskObj
        attenuationMaskObj
    end
    
    properties (Hidden)%, Access=protected) % dave jan 2025 - trying to optimise over angle, might regret this change later
        Defocus
        pupilMask
        fieldBFP
        chirpZTransform
    end
    properties (Dependent)
        oversampling (1,1)
        positionInPixelFromCenter (1,3)
        positionInPixelFromOrigin (1,3)
        positionInNanometerFromCenter (1,3)
    end
    properties (Dependent, Hidden)
        unitObjectSpace (1,1)
        unitKSpace (1,1)
    end

    methods
        %% Constructor
        function obj = PSF(par)
            if nargin > 0
                obj = setInputParameters('PSF', obj, par);
            end
            obj = obj.createImage();
        end

        %% Get methods for dependent properties
        function oversampling = get.oversampling(obj)
            oversampling = size(obj.pixelSensitivityMask, 1);
        end
        function positionInPixelFromCenter = get.positionInPixelFromCenter(obj)
            positionInPixelFromCenter = obj.position.inPixels(obj.pixelSize);
        end
        function positionInPixelFromOrigin = get.positionInPixelFromOrigin(obj)
            pos = obj.positionInPixelFromCenter;
            positionInPixelFromOrigin = [(obj.nPixels+1)/2 + pos(1), (obj.nPixels+1)/2 + pos(2)];
        end
        function positionInNanometerFromCenter = get.positionInNanometerFromCenter(obj)
            positionInNanometerFromCenter = obj.positionInPixelFromCenter .* obj.pixelSize.inMeter .* 1e9;
        end
        function unitObjectSpace = get.unitObjectSpace(obj)
            unitObjectSpace = obj.pixelSize.inMeter / obj.oversampling;
        end
        function unitKSpace = get.unitKSpace(obj)
            % Largest spatial frequency passed by objective lens
            maxSpatialFrequency = obj.objectiveNA / obj.wavelength.inMeter;
            maxAngularFrequency = 2*pi * maxSpatialFrequency;
            % Unit in pupil space (k-space)
            unitKSpace = 2 * maxAngularFrequency / obj.nDiscretizationBFP; % full range covers 2 times the maximum frequency
        end

        %% Functions for calculation of PSF

        function obj = createImage(obj)
            [obj.Defocus, obj.pupilMask, obj.chirpZTransform, obj.phaseMaskObj, obj.attenuationMaskObj] = obj.setup();
            
            bfp = BackFocalPlane(obj);
            % bfp = BackFocalPlane_gaussian(obj); % dave feb 2025 - use if just want gaussian
            obj.backFocalPlane = bfp;

            % Apply phase mask
            obj.fieldBFP.x = obj.phaseMaskObj.apply(bfp.electricField.x);
            obj.fieldBFP.y = obj.phaseMaskObj.apply(bfp.electricField.y);

            % Apply attenuation mask
            obj.fieldBFP.x = obj.attenuationMaskObj.apply(obj.fieldBFP.x);
            obj.fieldBFP.y = obj.attenuationMaskObj.apply(obj.fieldBFP.y);

            psf = zeros(obj.nPixels,obj.nPixels); 
            for k=1:size(obj.stageDrift.motion,1)
                % Apply aberrations
                aberrations = getAberrations(obj,k);
                aberratedFieldBFP = applyAberrations(obj, aberrations);
                
                % Get image from BFP field
                psf = psf + getIntensitiesCamera(obj, aberratedFieldBFP)./size(obj.stageDrift.motion,1);
            end

            psf = adjustExcitation(obj, psf);
            psf = applyShotNoise(obj, psf);
            psf = addBackgroundNoise(obj, psf);
            obj.image = psf;
        end

        function [Defocus, pupilMask, chirpZTransform, phaseMaskObj, attenuationObj] = setup(obj)
            Defocus = sphericalAberrationFromRefractiveIndexMismatch(obj);
            parMask.nGrid = obj.nDiscretizationBFP;
            pupilMask = Mask(parMask);
            chirpZTransform = ChirpZTransform(obj);
            phaseMaskObj = obj.phaseMask(obj.nDiscretizationBFP);
            attenuationObj = obj.attenuation(obj.nDiscretizationBFP);
        end

        function aberrations = getAberrations(obj,k)
            zernikeConstant = obj.unitObjectSpace * obj.unitKSpace * obj.nDiscretizationBFP /4;
            pos = obj.position.inPixels(obj.pixelSize);
            relativeMotion = obj.stageDrift.motion.inPixels(obj.pixelSize);
            
            if size(obj.stageDrift.motion,2)==3
                axialMotionInMeter = obj.stageDrift.motion.inMeter;
                axialMotionInMeter = axialMotionInMeter(:,3);
            else
                axialMotionInMeter = zeros(size(obj.stageDrift.motion,1),1); 
            end
           
            % (circle) Zernike polynomials
            Z = ZernikePolynomials.getInstance(obj.pupilMask);
            coeffPos = - (pos(1:2) + relativeMotion(k,1:2)) .* obj.oversampling .* zernikeConstant / (2*pi); % x/y-tilt
            if obj.astigmatism ~= 0
                coeffAstigmatism = obj.astigmatism; % vertical astigmatism
                aberrations = Z.getAberration([2,3,6],[coeffPos,coeffAstigmatism]) + 1/(2*pi)*(axialMotionInMeter(k)+obj.defocus.inMeter)*obj.Defocus;
            else
                aberrations = Z.getAberration([2,3],coeffPos) + 1/(2*pi)*(axialMotionInMeter(k)+obj.defocus.inMeter)*obj.Defocus;
            end
        end

        function fieldBFP = applyAberrations(obj, aberrations)
            mask = obj.pupilMask.values .* exp( 1i*2*pi*aberrations) ;%+ 1i*(obj.defocus.inMeter)*obj.Defocus );
            fieldBFP.x = obj.fieldBFP.x .* mask;
            fieldBFP.y = obj.fieldBFP.y .* mask;
        end

        function psf = adjustExcitation(obj, psf)
            % Dipole excitation probability dependent on
            % angle between electric field and dipole orientation
            if obj.reducedExcitation
                psf = psf * (sin(obj.dipole.inclination))^2;
            end
        end
        function psf = applyShotNoise(obj, psf)
            if obj.shotNoise
                psf = poissrnd(psf);
            end
        end
        function psf = addBackgroundNoise(obj, psf)
            if obj.backgroundNoise ~= 0
                psf = psf + poissrnd(obj.backgroundNoise*ones(size(psf)));
            end
        end

        % Additional functions defined in separate files
        [SA_out,Defocus,a] = sphericalAberrationFromRefractiveIndexMismatch(obj, removeDefocus)
        [psf, I_xx, I_yx] = getIntensitiesCamera(obj, mask)
        
        par = readParameters(obj)
        
        
        %% Plot function
        function plot(obj,setLimitsFullRange)
            if nargin==2 && setLimitsFullRange
                imagesc(obj.image(:,:),[0 max(obj.image(:))]);
            else
                imagesc(flip(obj.image(:,:)));
            end
            colorbar;
            axis equal; axis tight; axis xy;
            title(['Pixelsize = ', num2str(obj.pixelSize.inNanometer), 'nm. ', ...
                'Dipole = [', num2str(obj.dipole.inclination),',',num2str(obj.dipole.azimuth) ,']']);
        end
    end
end

function [Defocus, SA, a] = sphericalAberrationFromRefractiveIndexMismatch(obj, removeDefocus)

% Outputs spherical aberration arising from a RI-mismatch, for z=1m inside
% the material of RI2. Output can be scaled to the existing z-value
% returned defocus phase is defined for medium n2.
% According to paper: Jesacher et al. 2010. Opex
% The "defocus-term" is removed analytically.
%
% Inputs:
% n1 ... refractive index of immersion medium
% n2 ... refractive index of target material
% NA ... numerical aperture of objective
% lambda_0 ... vacuum wavelength
% a  ... coefficient of "Defocus" contained in SA
% rd ... 0 or 1, if set to 1, defocus is removed

N = obj.nDiscretizationBFP;
NA = obj.objectiveNA;
n2 = obj.refractiveIndices(3);
k = 2*pi / obj.wavelength.inMeter;

par.nGrid = N;
par.spacing = 2/(N-1);
par.mode = 'FFT';
mask = Mask(par);
[normR, ~] = getPolarCoordinates(mask);

% Mean value of spherical defocus
MW = @(RI) 2/3*k*(-(-1 + RI^2/NA^2)^(3/2)+(RI^2/NA^2)^(3/2))*NA;
% Defocus function
Def = @(RI) real(k*NA*sqrt(RI^2/NA^2-normR.^2)-MW(RI)).*mask.values; %spherical defocus in medium with refractive index RI (for 1m of defocus)
Defocus = Def(n2);

if nargout > 1
    n1 = obj.refractiveIndices(3);
    SA = -(Def(n2)-Def(n1)); %spherical aberration phase (complete, i.e including defocus) (for 1m of defocus), (see paper: Jesacher & Booth, Opex 2010)

    a_spher = -1/72*k^2*pi*(72*n2^2-36*NA^2+(32*(-(-1 + n2^2/NA^2)^(3/2) + n2^3/NA^3)*(n1^3-n1^2*sqrt((n1 - NA)*(n1 + NA))+NA^2*sqrt((n1 - NA)*(n1 + NA))))/NA - (32*(n2^3 - n2^2*sqrt((n2 - NA)*(n2 + NA))+NA^2*sqrt((n2 - NA)*(n2 + NA)))^2)/NA^4+(9/NA^2)*(2*(-n1^3*n2 - n1*n2^3 + n1^2*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))+sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))*(n2^2-2*NA^2)) - (n1^2 - n2^2)^2*(log((n1 - n2)^2)-log(n1^2 + n2^2 - 2*(NA^2 + sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))))));
    def_norm = -((k^2*(16*n2^6 - 24*n2^4*NA^2 + 6*n2^2*NA^4 + NA^6 - 16*n2^5*sqrt(n2^2 - NA^2) + 16*n2^3*NA^2*sqrt(n2^2 - NA^2))*pi)/(18*NA^4));
    a = a_spher/def_norm; %coefficient of Defocus contained in SA (without normalization)

    if removeDefocus
        SA = SA - a_spher*Def(n2)/def_norm;
    end
end

endclassdef ZernikePolynomials < handle

    properties (Access=private)
        polynomials % ordering = Noll indices
        mask Mask
    end

    methods (Static)
        function obj = getInstance(mask)
            persistent instance
            if isempty(instance) || ~isequal(size(instance.getMask), size(mask.values)) || ~all(instance.getMask == mask.values,'all')
                instance = ZernikePolynomials(mask);
            end
            obj = instance;
        end
    end

    methods (Access=private)
        function obj = ZernikePolynomials(mask)
            obj.mask = mask;
            [m,n] = size(mask.values);
            obj.polynomials = NaN(m,n,0);
        end
    end

    methods (Access=public)
        function polynomials = getPolynomials(obj)
            polynomials = obj.polynomials;
        end

        function mask = getMask(obj)
            mask = obj.mask.values;
        end

        function numberCalculatedPolynomials = getNumberCalculatedPolynomials(obj)
            numberCalculatedPolynomials = size(obj.polynomials, 3);
        end

        function aberration = getAberration(obj, indicesNoll, coefficients)
            assert(numel(indicesNoll) == numel(coefficients))

            maxIndexRequired = max(indicesNoll);
            numberAlreadyCalculated = obj.getNumberCalculatedPolynomials;
            if (maxIndexRequired > numberAlreadyCalculated)
                obj.polynomials = cat(3, obj.polynomials, ...
                    obj.calculatePolynomials((numberAlreadyCalculated+1):maxIndexRequired));
            end
            % Weighted sum of Zernike polynomials
            weightedPolynomials = multiplyPagewise(obj.polynomials(:,:,indicesNoll), coefficients);
            aberration = sum( weightedPolynomials,3 );
        end

        function polynomials = calculatePolynomials(obj, nollIndices)
            [n,m] = ZernikePolynomials.noll2ZernikeIndices(nollIndices);
            polynomials = ZernikePolynomials.calculateZernikePolynomials( n, m, obj.mask, 1 );
        end
    end


    methods (Static)
        function [n,m] = noll2ZernikeIndices(j)
            n = floor(sqrt(2.*j) - 1/2);
            s = mod(n,2);
            m = (1-2.*mod(j,2)) .* (2 .* floor((1+2.*j-n.*(n+1)+s) ./ 4) - s);
        end

        function polynomials = calculateZernikePolynomials(n, m, mask, doNormalize)
            %   The radial Zernike polynomials are the radial portion of the
            %   Zernike functions, which are an orthogonal basis on the unit
            %   circle.  The series representation of the radial Zernike
            %   polynomials is
            %
            %          (n-m)/2
            %            __
            %    m      \       s                                          n-2s
            %   Z(r) =  /__ (-1)  [(n-s)!/(s!((n-m)/2-s)!((n+m)/2-s)!)] * r
            %    n      s=0
            %
            %   The following table shows the first 12 polynomials.
            %
            %       n    m    Zernike polynomial    Normalization
            %       ---------------------------------------------
            %       0    0    1                        sqrt(2)
            %       1    1    r                        sqrt(4)
            %       2    0    2*r^2 - 1                sqrt(6)
            %       2    2    r^2                      sqrt(6)
            %       3    1    3*r^3 - 2*r              sqrt(8)
            %       3    3    r^3                      sqrt(8)
            %       4    0    6*r^4 - 6*r^2 + 1        sqrt(10)
            %       4    2    4*r^4 - 3*r^2            sqrt(10)
            %       4    4    r^4                      sqrt(10)
            %       5    1    10*r^5 - 12*r^3 + 3*r    sqrt(12)
            %       5    3    5*r^5 - 4*r^3            sqrt(12)
            %       5    5    r^5                      sqrt(12)
            %       ---------------------------------------------

            %% Check input
            if nargin<5
                doNormalize = 1;
            end

            assert( all(n>=abs(m)) && all(abs(m)>=0), 'Check input! It must hold that n>=|m|>=0.')
            %assert( all(0<=rho(:)) && all(rho(:)<=1), 'Radius must be between 0 and 1!')

            %% Vectorize
            sizeMask = size(mask.values);
            [rho,theta] = mask.getPolarCoordinates();
            rho = rho(:);
            theta = theta(:);
            mAbs = abs(m);

            %% Calculate polynomials

            [requiredPowersRho, powerIdx] = ZernikePolynomials.precalculateRequiredPowersOfRho(rho, n, m);

            polynomials = NaN(sizeMask(1),sizeMask(2),numel(n));
            for j = 1:numel(n)

                A = (n(j)+mAbs(j))/2;
                D = (n(j)-mAbs(j))/2;

                %% Radiual polynomials (sum from k = 0 to (n-m)/2)

                maxFactorial = max([n(j),A]);
                cumProds = cumprod([1,1:maxFactorial]);

                % (-1)^k (n-k)!
                nominator = (-1).^(0:D) .*  cumProds(n(j)-(0:D)+1);

                % k! * ((n+m)/2-k)! * ((n-m)/2-k)!
                F = cumProds((0:D)+1);
                denominator =  F .* cumProds(A-(0:D)+1) .* fliplr(F);

                % r^(n-2k)
                powers = n(j)-2*(0:D);
                powersRho = requiredPowersRho(:,powerIdx(powers+1));

                Rnm = sum(nominator./denominator .* powersRho,2);

                %% Zernike polynomials
                if m(j)==0
                    Z = Rnm;
                elseif m(j)>0 % 'even' Zernike polynomials
                    Z = Rnm .* cos(theta*mAbs(j));
                else % 'odd' Zernike polynomials
                    Z = Rnm .* sin(theta*mAbs(j));
                end

                %% Normalization
                if doNormalize
                    if m(j)==0
                        Z = Z *sqrt(n(j)+1);
                    else
                        Z = Z * sqrt(2*(n(j)+1));
                    end
                end

                polynomials(:,:,j) = reshape(Z,sizeMask).*mask.values;
            end
        end

        function [requiredPowersRho, powerIdx] = precalculateRequiredPowersOfRho(rho, n, m)
            isEven = mod(n,2);
            if all(isEven) || all(~isEven)
                requiredPowers = min(abs(m)):2:max(n);
            else
                requiredPowers = min(abs(m)):1:max(n);
            end
            if requiredPowers(1)==0 % faster version if power 0 is required
                requiredPowersRho = arrayfun(@(p) rho.^p, requiredPowers(2:end), 'UniformOutput', false);
                requiredPowersRho = cat(2, requiredPowersRho{:});
                requiredPowersRho = [ones(length(rho),1) requiredPowersRho];
            else
                requiredPowersRho = arrayfun(@(p) rho.^p, requiredPowers, 'UniformOutput',false);
                requiredPowersRho = cat(2, requiredPowersRho{:});
            end
            powerIdx = NaN(max(n)+1,1);
            powerIdx(requiredPowers+1) = 1:length(requiredPowers);
        end
    end
end
