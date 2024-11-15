% function patches = extract_patches(centroids_image_coords, current_frame_number, patch_width, image)
% % EXTRACT_PATCHES Extracts patches from the image based on centroids.
% % 
% % INPUTS:
% %   centroids_image_coords - Nx3 matrix with columns [frame, x, y]
% %   current_frame_number   - Frame number to process
% %   patch_width            - Width of the square patch (in pixels)
% %   image                  - Input image (2D matrix)
% %
% % OUTPUT:
% %   patches - Cell array containing extracted patches
% 
%     % Initialize an empty cell array to store patches
%     patches = {};
% 
%     % Loop over each centroid in centroids_image_coords
%     for k = 1:size(centroids_image_coords, 1)
%         f = centroids_image_coords(k, 1);
%         x = centroids_image_coords(k, 2);
%         y = centroids_image_coords(k, 3);
% 
%         % Only operate on the current frame number
%         if f == current_frame_number
%             % Define the coordinates for the patch
%             x_start = x - floor(patch_width / 2);
%             x_end = x + floor(patch_width / 2);
%             y_start = y - floor(patch_width / 2);
%             y_end = y + floor(patch_width / 2);
% 
%             % Determine the padding needed for each side
%             pad_left = max(0, 1 - x_start);
%             pad_right = max(0, x_end - size(image, 2));
%             pad_top = max(0, 1 - y_start);
%             pad_bottom = max(0, y_end - size(image, 1));
% 
%             % Apply padding to the image (using zero padding)
%             padded_image = padarray(image, [pad_top, pad_left], 0, 'pre');
%             padded_image = padarray(padded_image, [pad_bottom, pad_right], 0, 'post');
% 
%             % Recalculate the new patch coordinates based on the padded image
%             x_start_padded = x_start + pad_left;
%             x_end_padded = x_end + pad_left;
%             y_start_padded = y_start + pad_top;
%             y_end_padded = y_end + pad_top;
% 
%             % Extract the patch from the padded image
%             patch = padded_image(y_start_padded:y_end_padded, x_start_padded:x_end_padded);
% 
%             % Store the patch in the cell array
%             patches{end+1} = patch;
%         end
%     end
% end

function [patches, outlined_image] = extract_patches(centroids_image_coords, current_frame_number, patch_width, image)
% EXTRACT_PATCHES Extracts patches from the image based on centroids and shows outlines and centroids.
% 
% INPUTS:
%   centroids_image_coords - Nx3 matrix with columns [frame, x, y]
%   current_frame_number   - Frame number to process
%   patch_width            - Width of the square patch (in pixels)
%   image                  - Input image (2D matrix)
%
% OUTPUTS:
%   patches        - Cell array containing extracted patches
%   outlined_image - The original image with red outlines and green centroids

    % Initialize an empty cell array to store patches
    patches = {};
    
    % Convert image to RGB if it is grayscale
    if size(image, 3) == 1
        outlined_image = repmat(image, [1, 1, 3]);
    else
        outlined_image = image;
    end

    % Loop over each centroid in centroids_image_coords
    for k = 1:size(centroids_image_coords, 1)
        f = centroids_image_coords(k, 1);
        x = centroids_image_coords(k, 2);
        y = centroids_image_coords(k, 3);

        % Only operate on the current frame number
        if f == current_frame_number
            % Define the coordinates for the patch
            x_start = x - round(patch_width / 2);
            y_start = y - round(patch_width / 2);
            width = patch_width;
            height = patch_width;

            % Ensure coordinates are within bounds
            if x_start < 1, x_start = 1; end
            if y_start < 1, y_start = 1; end

            % Extract the patch
            x_end = min(x_start + width - 1, size(image, 2));
            y_end = min(y_start + height - 1, size(image, 1));
            patch = image(y_start:y_end, x_start:x_end);

            % Store the patch in the cell array
            patches{end+1} = patch;

            % Draw a red rectangle on the original image for the patch
            outlined_image = insertShape(outlined_image, 'Rectangle', ...
                                         [x_start, y_start, width, height], ...
                                         'Color', 'red', 'LineWidth', 2);
            
        end
    end
end
