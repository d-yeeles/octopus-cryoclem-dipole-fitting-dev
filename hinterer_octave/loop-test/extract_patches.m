function [patches, patch_centres_x_nm, patch_centres_y_nm, outlined_image] = extract_patches( ...
    frameNumber_array, ...
    centroidsX_array, ...
    centroidsY_array, ...
    current_frame_number, ...
    pixel_width_nm, ...
    patch_width_nm, ...
    image_width_nm, ...
    image_width_px, ...
    image ...
    )

% EXTRACT_PATCHES Extracts patches from the image based on centroids and shows outlines and centroids.
%
% OUTPUTS:
%   patches        - cell array containing extracted patches
%   patch_centres  - cell array containing patch centres in nm, which may differ
%   from the input coord arrays because of cases when patch would go over
%   the edge of the image
%   outlined_image - original image with patches drawn on

    % Initialize an empty cell array to store patches
    patches = {};
    patch_centres_x_nm = [];
    patch_centres_y_nm = [];
    
    % Convert image to RGB if it is grayscale
    if size(image, 3) == 1
        outlined_image = repmat(image, [1, 1, 3]);
    else
        outlined_image = image;
    end

    % Loop over each centroid in centroids_image_coords
    for i = 1:length(centroidsX_array)
        f = frameNumber_array(i);
        x_nm = centroidsX_array(i);
        y_nm = centroidsY_array(i);

        % Only operate on the current frame number
        if f == current_frame_number

            % Nominal patch boundaries, assuming no going over the edges
            x_start_nm = x_nm - patch_width_nm/2;
            y_start_nm = y_nm - patch_width_nm/2;
            x_end_nm = x_nm + patch_width_nm/2;
            y_end_nm = y_nm + patch_width_nm/2;



            % Handle patches that would go over the edge

            % if left/top would go outside image, set it to edge instead
            if x_start_nm < -image_width_nm/2
                x_start_nm = -image_width_nm/2;
                x_end_nm = x_start_nm + patch_width_nm; % the -1 to account for matlab 1 indexing
            end
            if y_start_nm < -image_width_nm/2
                y_start_nm = -image_width_nm/2;
                y_end_nm = y_start_nm + patch_width_nm;
            end

            % if right/bottom would go outside image, set it to image width/height instead
            % and retroactively make it patch width wide/tall
            if x_end_nm > image_width_nm/2
                x_end_nm = image_width_nm/2;
                x_start_nm = x_end_nm - patch_width_nm; % the +1 to account for matlab 1 indexing
            end
            if y_end_nm > image_width_nm/2
                y_end_nm = image_width_nm/2;
                y_start_nm = y_end_nm - patch_width_nm; % the +1 to account for matlab 1 indexing
            end



            patch_centre_x_nm = x_start_nm + (x_end_nm - x_start_nm)/2;
            patch_centre_y_nm = y_start_nm + (y_end_nm - y_start_nm)/2;

            % convert things to px to define patch on image
            patch_width_px = roundToOdd(patch_width_nm/pixel_width_nm);

            % patch_centre_x_px = nm_to_px(patch_centre_x, image_width_nm, image_width_px, 'x');
            % patch_centre_y_px = nm_to_px(patch_centre_y, image_width_nm, image_width_px, 'y');
            % x_px = nm_to_px(x, image_width_nm, image_width_px, 'x');
            % y_px = nm_to_px(y, image_width_nm, image_width_px, 'y');

            % patch_centre_x_px = patch_centre_x_nm/51.2 + (image_width_px + 1)/2;
            % patch_centre_y_px = image_width_px - (patch_centre_y_nm/51.2 + (image_width_px + 0)/2);
            % x_px = x_nm/51.2 + (image_width_px + 1)/2;
            % y_px = image_width_px - (y_nm/51.2 + (image_width_px + 0)/2);
            patch_centre_x_px = nm_to_px(patch_centre_x_nm, pixel_width_nm, image_width_px, 'x');
            patch_centre_y_px = nm_to_px(patch_centre_y_nm, pixel_width_nm, image_width_px, 'y');
            x_px = nm_to_px(x_nm, pixel_width_nm, image_width_px, 'x');
            y_px = nm_to_px(y_nm, pixel_width_nm, image_width_px, 'y');

            x_start_px = patch_centre_x_px - (patch_width_px - 1)/2;
            y_start_px = patch_centre_y_px - (patch_width_px - 1)/2;
            x_end_px = patch_centre_x_px + (patch_width_px - 1)/2;
            y_end_px = patch_centre_y_px + (patch_width_px - 1)/2;

            patch = image(int32(round(y_start_px)):int32(round(y_end_px)), int32(round(x_start_px)):int32(round(x_end_px)));

            % Store the patch in the cell array
            patch_centres_x_nm(end+1) = patch_centre_x_nm;
            patch_centres_y_nm(end+1) = patch_centre_y_nm;
            patches{end+1} = patch;

            red_value = 1 - (i - 1) / (length(centroidsX_array) - 1); % Red decreases
            blue_value = (i - 1) / (length(centroidsX_array) - 1);    % Blue increase            
            current_color = uint8([red_value, 0, blue_value] * 255);
            if length(centroidsX_array) == 1
                current_color = 'red';
            end

            % Draw a red rectangle on the original image for the patch
            outlined_image = insertShape(outlined_image, 'Rectangle', ...
                                         [x_start_px, y_start_px, patch_width_px, patch_width_px], ...
                                         'Color', current_color, 'LineWidth', 1);
            outlined_image = insertShape(outlined_image, 'Circle', [x_px, y_px, 0], 'Color', 'red', 'LineWidth', 1);
            
        end
    end
end
