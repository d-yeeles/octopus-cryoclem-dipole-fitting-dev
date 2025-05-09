function run_fit_stack_commandline_thunderstorm_simulated(varargin)
    if nargin < 4
        error('Usage: matlab -nodisplay -r "run_fit_stack_commandline_thunderstorm_simulated(''input_image_and_thunderstorm_path'', ''model'', ''patch_width_nm'', ''ROI_centre_x_nm'', ''ROI_centre_y_nm'', ''ROI_width_nm'', ''ROI_height_nm'', ''starting_frame_index'', ''ending_frame_index'')"');
    end
    
    input_image_and_thunderstorm_path = strip_quotes(varargin{1});
    model = strip_quotes(varargin{2});
    patch_width_nm_str = strip_quotes(varargin{3});
    ROI_centre_x_nm_str = strip_quotes(varargin{4});
    ROI_centre_y_nm_str = strip_quotes(varargin{5});
    ROI_width_nm_str = strip_quotes(varargin{6});
    ROI_height_nm_str = strip_quotes(varargin{7});
    starting_frame_index = strip_quotes(varargin{8});
    ending_frame_index = strip_quotes(varargin{9});

    ROI_centre_x_nm = str2double(ROI_centre_x_nm_str);
    ROI_centre_y_nm = str2double(ROI_centre_y_nm_str);
    ROI_width_nm = str2double(ROI_width_nm_str);
    ROI_height_nm = str2double(ROI_height_nm_str);
    patch_width_nm = str2double(patch_width_nm_str);
    starting_frame_index = str2double(starting_frame_index);
    ending_frame_index = str2double(ending_frame_index);

    fit_stack_commandline_thunderstorm_simulated(input_image_and_thunderstorm_path, model, patch_width_nm, ROI_centre_x_nm, ROI_centre_y_nm, ROI_width_nm, ROI_height_nm, starting_frame_index, ending_frame_index);
    exit;
end

function cleanStr = strip_quotes(str)
    if ischar(str) || isstring(str)
        % remove any surrounding quotes
        if (length(str) >= 2)
            if (str(1) == '"' && str(end) == '"') || (str(1) == '''' && str(end) == '''')
                str = str(2:end-1);
            end
        end
    end
    cleanStr = str;
end