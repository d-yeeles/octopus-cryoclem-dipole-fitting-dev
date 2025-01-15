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
        
            % Extracting parameters from psfObj
            RI = psfObj.refractiveIndices;
            dipole = psfObj.dipole;
            hIntermediate = psfObj.heightIntermediateLayer.inMeter;
            pos = psfObj.position.inMeter;
            z = pos(3);
            focalLength = psfObj.objectiveFocalLength.inMeter;
        
            % Coordinates in the objective pupil
            par.nGrid = obj.nGrid;
            par.spacing = obj.unitKSpace;
            pupilMask = Mask(par);
            [Kr, ~] = pupilMask.getPolarCoordinates; % Radial coordinates in the pupil
        
            % Define Gaussian parameters
            sigma = 0.5; % Standard deviation of the Gaussian
            amplitude = 1; % Peak amplitude of the Gaussian
        
            % Gaussian Electric Field
            GaussianField = amplitude * exp(-Kr.^2 / (2 * sigma^2));
        
            % Output Gaussian electric field for both x and y components
            E_BFP.x = GaussianField; % x-component of the electric field
            E_BFP.y = GaussianField; % y-component of the electric field (same as x for symmetry)
        
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