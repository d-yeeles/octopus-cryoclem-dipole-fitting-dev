function run_fit_stack_commandline_thunderstorm(varargin)
    if nargin < 5
        error('Usage: matlab -nodisplay -r "run_fit_stack_commandline_thunderstorm(''input_image_path'', ''input_thunderstorm_path'', ''output_results_path'', ''model'', ''patch_width_nm'', ''first_frame'', ''last_frame'')"');
    end
    
    input_image_path = strip_quotes(varargin{1});
    input_thunderstorm_path = strip_quotes(varargin{2});
    output_results_path = strip_quotes(varargin{3});
    model = strip_quotes(varargin{4});
    patch_width_nm_str = strip_quotes(varargin{5});
    first_frame_str = strip_quotes(varargin{6});
    last_frame_str = strip_quotes(varargin{7});

    patch_width_nm = str2double(patch_width_nm_str);
    first_frame = str2double(first_frame_str);
    last_frame = str2double(last_frame_str);

    % % If frame_array is a string of comma-separated numbers like "1,2,3,4"
    % if ischar(frame_array) || isstring(frame_array)
    %     % Remove any brackets if present
    %     frame_array = strrep(frame_array, '[', '');
    %     frame_array = strrep(frame_array, ']', '');
    %     % Split by comma
    %     frame_array_str = strsplit(frame_array, ',');
    %     % Convert each element to number
    %     frame_array = arrayfun(@str2double, frame_array_str);
    % end

    fit_stack_commandline_thunderstorm(input_image_path, input_thunderstorm_path, output_results_path, model, patch_width_nm, first_frame, last_frame);
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