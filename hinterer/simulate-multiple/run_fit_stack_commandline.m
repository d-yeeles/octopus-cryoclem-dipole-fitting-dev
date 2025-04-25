function run_fit_stack_commandline(varargin)

    if nargin < 5
        error('Usage: matlab -nodisplay -r "run_fit_multiple_psfs(''stack_dir'', ''stack_name'', ''results_path'', ''model'', ''patch_width_nm'')"');
    end
    
    stack_dir = varargin{1};
    stack_name = varargin{2};
    results_path = varargin{3};
    model = varargin{4};
    patch_width_nm = str2double(varargin{5});

    fit_stack_commandline(stack_dir, stack_name, results_path, model, patch_width_nm);
    
    exit;
end