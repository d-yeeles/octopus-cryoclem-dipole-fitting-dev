function nm = px_to_nm(px, pixel_size_nm, image_size_px, axis)
    % px_to_nm - Converts image coordinates (px) to physical coordinates (nm)
    %
    % Inputs:
    %   px              - Position in pixels from top-left (image coordinates)
    %   pixel_size_nm   - Size of one pixel in nanometers
    %   image_size_px   - Size of the image in pixels [width, height] or scalar if square
    %   axis            - 'x' or 'y' axis to specify coordinate conversion direction
    %
    % Outputs:
    %   nm              - Position in nanometers from image center (physical coordinates)
    %
    % Physical coordinates: Origin at center of image, positive right/up
    % Image coordinates: Origin at top-left, positive right/down
    
    % Handle scalar or vector input for image_size_px
    if length(image_size_px) == 1
        % Square image
        width_px = image_size_px;
        height_px = image_size_px;
    else
        % Rectangular image
        width_px = image_size_px(1);
        height_px = image_size_px(2);
    end
    
    % Get center of image in pixel coordinates
    center_x_px = floor(width_px/2) + 1;
    center_y_px = floor(height_px/2) + 1;
    
    % Convert based on axis
    if strcmpi(axis, 'x')
        % For X-axis: px → nm, shift origin and scale
        nm = (px - center_x_px) * pixel_size_nm;
    elseif strcmpi(axis, 'y')
        % For Y-axis: px → nm, shift origin and scale, but also flip direction
        % Decreasing px (up) → positive nm (up)
        nm = (center_y_px - px) * pixel_size_nm;
    else
        error('Axis must be either ''x'' or ''y''');
    end
end