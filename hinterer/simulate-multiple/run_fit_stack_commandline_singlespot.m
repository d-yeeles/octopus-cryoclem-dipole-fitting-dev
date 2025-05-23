function run_fit_stack_commandline_singlespot(varargin)
    if nargin < 5
        error('Usage: matlab -nodisplay -r "run_fit_stack_commandline_singlespot(''stack_path'', ''params_path'', ''results_path'', ''model'', ''patch_width_nm'', ''starting_frame_index'', ''ending_frame_index'')"');
    end
    
    stack_path = strip_quotes(varargin{1});
    params_path = strip_quotes(varargin{2});
    results_path = strip_quotes(varargin{3});
    model = strip_quotes(varargin{4});
    patch_width_nm_str = strip_quotes(varargin{5});
    starting_frame_index = strip_quotes(varargin{6});
    ending_frame_index = strip_quotes(varargin{7});

    patch_width_nm = str2double(patch_width_nm_str);
    starting_frame_index = str2double(starting_frame_index);
    ending_frame_index = str2double(ending_frame_index);

    fit_stack_commandline_singlespot(stack_path, params_path, results_path, model, patch_width_nm, starting_frame_index, ending_frame_index);
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