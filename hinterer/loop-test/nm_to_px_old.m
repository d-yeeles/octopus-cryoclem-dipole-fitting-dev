
% function x_px = nm_to_px(x_nm, image_width_nm, image_width_pixels, x_or_y_flag)
% 
%     if x_or_y_flag == 'x'
%         shifted_x = x_nm + image_width_nm / 2;  % Shift to start from the left
%         relative_x = shifted_x / image_width_nm; % Normalize to [0, 1]
%         x_px = relative_x * (image_width_pixels - 1) + 1;  % Convert to pixel coordinates and adjust for 1-based indexing. The 1.5 is so we map centre to centre, not centre to edge of pixel
%     end
% 
%     if x_or_y_flag == 'y'
%         shifted_x = x_nm + image_width_nm / 2;  % Shift to start from the top (center at (0,0))
%         relative_x = shifted_x / image_width_nm; % Normalize to [0, 1]
%         x_px = relative_x * (image_width_pixels - 1) + 1;  % Convert to pixel coordinates and adjust for 1-based indexing
%         x_px = image_width_pixels - x_px + 1;  % Invert y-axis: top to bottom flip
%     end
% 
% end


function coord_px = nm_to_px(coord_nm, pixel_width_nm, image_width_px, x_or_y_flag)

    if x_or_y_flag == 'x'
        coord_px = (coord_nm / pixel_width_nm) + (image_width_px + 1) / 2;
    end

    if x_or_y_flag == 'y'
        coord_px = image_width_px - (coord_nm / pixel_width_nm) - (image_width_px + 1) / 2 + 1;
    end

end
