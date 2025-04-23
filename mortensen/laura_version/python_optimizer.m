
function [x, fval] = python_optimizer(x0, lb, ub, options_json)
    % This function serves as a bridge between MATLAB's fmincon and Python
    % It will read/write to temporary files to communicate with Python
    
    % Parse options
    options = optimoptions(@fmincon, ...
      'Display', 'off', ... % Do not display output
      'StepTolerance', 1e-6, ...          % Stop when the step size is less than 1e-6
      'OptimalityTolerance', 1e-6, ...    % Stop when the gradient is less than 1e-6
      'FunctionTolerance', 1e-6 ...        % Stop if the function value change is less than 1e-6
    );
    
    % Create constraint function for the unit sphere constraint
    function [c, ceq] = sphere_constraint(x)
        c = [];
        ceq = x(1)^2 + x(2)^2 + x(3)^2 - 1;
    end
    
    % Create objective function that calls Python
    function f = objective(x)
        % Write parameters to temp file
        param_file = 'matlab_params_temp.json';
        fid = fopen(param_file, 'w');
        fprintf(fid, jsonencode(x'));
        fclose(fid);
        
        % Signal Python to evaluate
        ready_file = 'matlab_ready.txt';
        fid = fopen(ready_file, 'w');
        fprintf(fid, '1');
        fclose(fid);
        
        % Wait for Python to calculate and provide the value
        result_file = 'python_result_temp.json';
        while ~exist(result_file, 'file')
            pause(0.01);
        end
        
        % Read result
        fid = fopen(result_file, 'r');
        result_str = fscanf(fid, '%s');
        fclose(fid);
        
        % Parse result
        f = str2double(result_str);
        
        % Clean up
        if exist(result_file, 'file'), delete(result_file); end
    end
    
    % Run the optimization
    [x, fval] = fmincon(@objective, x0, [], [], [], [], lb, ub, @sphere_constraint, options);
end
        