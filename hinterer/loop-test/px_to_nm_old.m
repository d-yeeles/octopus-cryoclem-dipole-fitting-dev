% function x_nm = px_to_nm(x_px, image_width_nm, image_width_pixels)
%     % Ensure MATLAB's 1-based indexing is handled correctly
%     % Convert from 1-based pixel index to 0-based index for calculation
%     x_px_shifted = x_px - 1;  
% 
%     % Convert from pixel coordinates to relative position in [0,1]
%     relative_x = x_px_shifted / (image_width_pixels - 1);  % Normalize to [0, 1]
% 
%     % Convert from relative position to physical coordinates in nm
%     shifted_x = relative_x * image_width_nm;
% 
%     % Shift back to the physical coordinate system where (0,0) is the center
%     x_nm = shifted_x - image_width_nm / 2;
% end

% function x_nm = px_to_nm(x_px, image_width_nm, image_width_pixels, x_or_y_flag)
% 
%     if x_or_y_flag == 'x'
%         x_px_shifted = x_px - 1;  % Convert from 1-based to 0-based indexing
%         relative_x = x_px_shifted / (image_width_pixels - 1);  % Normalize to [0, 1]
%         shifted_x = relative_x * image_width_nm;  % Convert to physical coordinates
%         x_nm = shifted_x - image_width_nm / 2;  % Shift back to physical coordinates (centered at 0)
%     end
% 
%     if x_or_y_flag == 'y'
%         x_px_shifted = x_px - 1;  % Convert from 1-based to 0-based indexing
%         relative_x = x_px_shifted / (image_width_pixels - 1);  % Normalize to [0, 1]
%         shifted_x = relative_x * image_width_nm;  % Convert to physical coordinates
%         x_nm = shifted_x - image_width_nm / 2;  % Shift back to physical coordinates (centered at 0)
%         x_nm = -x_nm;  % Flip y-axis: invert the y-coordinate to match the physical system's positive-up orientation
%     end
% 
% end

function coord_nm = px_to_nm(coord_px, pixel_width_nm, image_width_px, x_or_y_flag)

    if x_or_y_flag == 'x'
        coord_nm = (coord_px - (image_width_px + 1) / 2) * pixel_width_nm;
    elseif x_or_y_flag == 'y'
        coord_nm = (image_width_px - coord_px - (image_width_px + 1) / 2 + 1) * pixel_width_nm;
    else
        error('Invalid x_or_y_flag. It must be either ''x'' or ''y''.');
    end

end

