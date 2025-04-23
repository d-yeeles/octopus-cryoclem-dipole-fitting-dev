% Coordinate Transformation Validation Script
% This script tests the reversibility of coordinate transformations
% between pixel and nanometer coordinate systems for different image sizes

% Clear workspace
clear all;
close all;

% Add path to make sure we can access the conversion functions
% addpath(genpath('../')); % Adjust if needed

% Define test parameters
pixel_size_nm = 100; % Size of a pixel in nm
test_positions_px = [10.5, 20.5; 50.25, 60.75; 100.33, 200.67]; % Test positions in pixels
image_sizes_px = [19, 115]; % Two different image sizes to test

% Print header
fprintf('Coordinate Transformation Validation\n');
fprintf('====================================\n\n');

% Test for each image size
for i = 1:length(image_sizes_px)
    image_size_px = image_sizes_px(i);
    
    fprintf('Testing with image size: %d x %d pixels\n', image_size_px, image_size_px);
    fprintf('----------------------------------------\n');
    
    % Get center of image in pixel coordinates (for reference)
    center_x_px = (image_size_px + 1) / 2;
    center_y_px = (image_size_px + 1) / 2;
    fprintf('Image center is at pixel coordinates: (%.2f, %.2f)\n\n', center_x_px, center_y_px);
    
    % Test each position
    for j = 1:size(test_positions_px, 1)
        original_x_px = test_positions_px(j, 1);
        original_y_px = test_positions_px(j, 2);
        
        % Calculate nm coordinates
        x_nm = px_to_nm(original_x_px, pixel_size_nm, image_size_px, 'x');
        y_nm = px_to_nm(original_y_px, pixel_size_nm, image_size_px, 'y');
        
        % Convert back to pixels
        recovered_x_px = nm_to_px(x_nm, pixel_size_nm, image_size_px, 'x');
        recovered_y_px = nm_to_px(y_nm, pixel_size_nm, image_size_px, 'y');
        
        % Calculate errors
        x_error_px = original_x_px - recovered_x_px;
        y_error_px = original_y_px - recovered_y_px;
        
        % Print results
        fprintf('Test position %d: (%.2f, %.2f) px\n', j, original_x_px, original_y_px);
        fprintf('  → Converted to nm: (%.2f, %.2f) nm\n', x_nm, y_nm);
        fprintf('  → Converted back to px: (%.2f, %.2f) px\n', recovered_x_px, recovered_y_px);
        fprintf('  → Error: (%.6f, %.6f) px\n\n', x_error_px, y_error_px);
    end
    
    % Now test the actual patch extraction and recentering logic
    fprintf('Testing patch extraction and recentering logic:\n');
    fprintf('--------------------------------------------\n');
    
    % Parameters similar to your original code
    patch_width_nm = 1000;
    patch_width_px = patch_width_nm/pixel_size_nm;
    patch_width_px = 2*floor(patch_width_px/2) + 1; % Enforce odd number
    half_width = floor(patch_width_px / 2);
    
    % Create a test position (simulated emitter)
    test_emitter_nm_x = 500; % nm from center
    test_emitter_nm_y = -300; % nm from center
    
    % Convert to pixel coordinates
    test_emitter_px_x = nm_to_px(test_emitter_nm_x, pixel_size_nm, image_size_px, 'x');
    test_emitter_px_y = nm_to_px(test_emitter_nm_y, pixel_size_nm, image_size_px, 'y');
    
    fprintf('Test emitter position: (%.2f, %.2f) nm → (%.2f, %.2f) px\n', test_emitter_nm_x, test_emitter_nm_y, test_emitter_px_x, test_emitter_px_y);
    
    % Simulate patch extraction logic from your code
    adjusted_centre_x_px = min(max(half_width + 1, round(test_emitter_px_x)), image_size_px - half_width);
    adjusted_centre_y_px = min(max(half_width + 1, round(test_emitter_px_y)), image_size_px - half_width);
    
    patch_start_x_px = adjusted_centre_x_px - half_width;
    patch_start_y_px = adjusted_centre_y_px - half_width;
    patch_end_x_px = patch_start_x_px + patch_width_px - 1;
    patch_end_y_px = patch_start_y_px + patch_width_px - 1;
    
    % Calculate actual patch center
    actual_patch_center_x_px = patch_start_x_px + (patch_width_px - 1)/2;
    actual_patch_center_y_px = patch_start_y_px + (patch_width_px - 1)/2;
    
    % Convert back to nm
    actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_size_px, 'x');
    actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_size_px, 'y');
    
    % Calculate expected offset (position relative to patch center)
    expected_offset_x_nm = test_emitter_nm_x - actual_patch_center_x_nm;
    expected_offset_y_nm = test_emitter_nm_y - actual_patch_center_y_nm;
    
    % Now reverse the process (simulate what happens after fitting)
    fit_result_x_nm = expected_offset_x_nm; % Assume perfect fit for testing
    fit_result_y_nm = expected_offset_y_nm;
    
    final_position_x_nm = fit_result_x_nm + actual_patch_center_x_nm;
    final_position_y_nm = fit_result_y_nm + actual_patch_center_y_nm;
    
    % Calculate error
    position_error_x_nm = test_emitter_nm_x - final_position_x_nm;
    position_error_y_nm = test_emitter_nm_y - final_position_y_nm;
    
    fprintf('Adjusted patch center: (%.2f, %.2f) px → (%.2f, %.2f) nm\n', actual_patch_center_x_px, actual_patch_center_y_px, actual_patch_center_x_nm, actual_patch_center_y_nm);
    
    fprintf('Expected offset within patch: (%.2f, %.2f) nm\n', expected_offset_x_nm, expected_offset_y_nm);
    
    fprintf('Final reconstructed position: (%.2f, %.2f) nm\n', final_position_x_nm, final_position_y_nm);
    
    fprintf('Position error: (%.6f, %.6f) nm\n\n', position_error_x_nm, position_error_y_nm);
end

% Additional test: check behavior when emitter is near image edge
fprintf('Edge Case Test: Emitter Near Image Edge\n');
fprintf('======================================\n');

for i = 1:length(image_sizes_px)
    image_size_px = image_sizes_px(i);
    edge_margin = 3; % pixels from edge
    
    % Create test positions near edges
    edge_positions = [
        edge_margin, edge_margin;                           % Top-left
        edge_margin, image_size_px - edge_margin;           % Bottom-left
        image_size_px - edge_margin, edge_margin;           % Top-right
        image_size_px - edge_margin, image_size_px - edge_margin; % Bottom-right
    ];
    
    fprintf('\nTesting with image size: %d x %d pixels\n', image_size_px, image_size_px);
    
    for j = 1:size(edge_positions, 1)
        edge_x_px = edge_positions(j, 1);
        edge_y_px = edge_positions(j, 2);
        
        % Convert to nm
        edge_x_nm = px_to_nm(edge_x_px, pixel_size_nm, image_size_px, 'x');
        edge_y_nm = px_to_nm(edge_y_px, pixel_size_nm, image_size_px, 'y');
        
        % Simulate patch extraction
        half_width = floor(patch_width_px / 2);
        adjusted_centre_x_px = min(max(half_width + 1, round(edge_x_px)), image_size_px - half_width);
        adjusted_centre_y_px = min(max(half_width + 1, round(edge_y_px)), image_size_px - half_width);
        
        % Check if adjustment happened
        x_adjusted = (adjusted_centre_x_px ~= round(edge_x_px));
        y_adjusted = (adjusted_centre_y_px ~= round(edge_y_px));
        
        fprintf('Edge position %d: (%.2f, %.2f) px → (%.2f, %.2f) nm\n', j, edge_x_px, edge_y_px, edge_x_nm, edge_y_nm);
        
        fprintf('  Adjusted center: (%.2f, %.2f) px [%s, %s]\n', adjusted_centre_x_px, adjusted_centre_y_px, iif(x_adjusted, 'ADJUSTED', 'unchanged'), iif(y_adjusted, 'ADJUSTED', 'unchanged'));
        
        % Calculate actual patch boundaries
        patch_start_x_px = adjusted_centre_x_px - half_width;
        patch_start_y_px = adjusted_centre_y_px - half_width;
        patch_end_x_px = patch_start_x_px + patch_width_px - 1;
        patch_end_y_px = patch_start_y_px + patch_width_px - 1;
        
        fprintf('  Patch boundaries: (%.2f, %.2f) to (%.2f, %.2f) px\n', patch_start_x_px, patch_start_y_px, patch_end_x_px, patch_end_y_px);
        
        % Calculate actual patch center
        actual_patch_center_x_px = patch_start_x_px + (patch_width_px - 1)/2;
        actual_patch_center_y_px = patch_start_y_px + (patch_width_px - 1)/2;
        actual_patch_center_x_nm = px_to_nm(actual_patch_center_x_px, pixel_size_nm, image_size_px, 'x');
        actual_patch_center_y_nm = px_to_nm(actual_patch_center_y_px, pixel_size_nm, image_size_px, 'y');
        
        % Calculate expected offset
        expected_offset_x_nm = edge_x_nm - actual_patch_center_x_nm;
        expected_offset_y_nm = edge_y_nm - actual_patch_center_y_nm;
        
        % Simulate reconstruction
        final_position_x_nm = expected_offset_x_nm + actual_patch_center_x_nm;
        final_position_y_nm = expected_offset_y_nm + actual_patch_center_y_nm;
        
        % Calculate error
        position_error_x_nm = edge_x_nm - final_position_x_nm;
        position_error_y_nm = edge_y_nm - final_position_y_nm;
        
        fprintf('  Position error: (%.6f, %.6f) nm\n\n', position_error_x_nm, position_error_y_nm);
    end
end

% Helper function for conditional formatting
function result = iif(condition, true_value, false_value)
    if condition
        result = true_value;
    else
        result = false_value;
    end
end