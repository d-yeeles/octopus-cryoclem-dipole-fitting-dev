function px = nm_to_px(nm, pixel_size_nm, image_size_px, axis)
    % nm_to_px - Converts physical coordinates (nm) to image coordinates (px)
    %
    % Inputs:
    %   nm              - Position in nanometers from image center (physical coordinates)
    %   pixel_size_nm   - Size of one pixel in nanometers
    %   image_size_px   - Size of the image in pixels [width, height] or scalar if square
    %   axis            - 'x' or 'y' axis to specify coordinate conversion direction
    %
    % Outputs:
    %   px              - Position in pixels from top-left (image coordinates)
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
        % For X-axis: nm → px, just shift origin and scale
        % Positive nm (right) → increasing px
        px = (nm / pixel_size_nm) + center_x_px;
    elseif strcmpi(axis, 'y')
        % For Y-axis: nm → px, shift origin and scale, but also flip direction
        % Positive nm (up) → decreasing px (down)
        px = center_y_px - (nm / pixel_size_nm);
    else
        error('Axis must be either ''x'' or ''y''');
    end
end