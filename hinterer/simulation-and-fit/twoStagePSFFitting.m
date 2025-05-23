function fitResult = twoStagePSFFitting(psfInit, parEst, model, Nlive_simple, Nlive_complex, tolerance_simple, tolerance_complex)
    % Two-stage adaptive PSF fitting with customizable live points and tolerances
    % 
    % INPUTS:
    %   psfInit, parEst, model - standard PSF fitting inputs
    %   Nlive_simple - live points for simple cases (optional, default: 200)
    %   Nlive_complex - live points for complex cases (optional, default: 400)
    %   tolerance_simple - tolerance for simple cases (optional, default: 0.08)
    %   tolerance_complex - tolerance for complex cases (optional, default: 0.03)
    
    % Set defaults if not provided
    if nargin < 4 || isempty(Nlive_simple)
        Nlive_simple = 200;
    end
    if nargin < 5 || isempty(Nlive_complex)
        Nlive_complex = 400;
    end
    if nargin < 6 || isempty(tolerance_simple)
        tolerance_simple = 0.08;
    end
    if nargin < 7 || isempty(tolerance_complex)
        tolerance_complex = 0.03;
    end
    
    fprintf('=== Two-Stage PSF Fitting ===\n');
    fprintf('Settings: Simple (Nlive=%d, tol=%.3f), Complex (Nlive=%d, tol=%.3f)\n', ...
        Nlive_simple, tolerance_simple, Nlive_complex, tolerance_complex);
    
    % First try: Always use simple settings
    fprintf('Initial Fit: Using standard sampling (Nlive=%d, tol=%.3f)\n', Nlive_simple, tolerance_simple);
    tic;
    
    % Run fit with simple settings
    fitResult = FitPSF_posterior_twostage(psfInit, parEst, model, false);
    fitResult = fitResult.runFittingWithParams(...
        'Nlive', Nlive_simple, ...
        'tolerance', tolerance_simple, ...
        'plotPosteriors', false);
    
    simple_time = toc;
    fprintf('  Initial fit completed in %.1f seconds\n', simple_time);
    
    % Check complexity of the result
    is_complex = analyzeComplexity(fitResult, psfInit, parEst);
    
    % Create a settings structure we'll add to the object's nestedSamplingParams
    timing_info = struct();
    timing_info.simple_time = simple_time;
    timing_info.total_time = simple_time;
    timing_info.complexity_detected = is_complex;
    
    % If complex, re-run with more live points
    if is_complex
        fprintf('Re-running with robust sampling (Nlive=%d, tol=%.3f) for complex posterior\n', ...
            Nlive_complex, tolerance_complex);
        tic;
        
        % Run final fit with complex settings
        fitResult = FitPSF_posterior_twostage(psfInit, parEst, model, false);
        fitResult = fitResult.runFittingWithParams(...
            'Nlive', Nlive_complex, ...
            'tolerance', tolerance_complex, ...
            'plotPosteriors', false);
        
        complex_time = toc;
        fprintf('  Complex fit completed in %.1f seconds\n', complex_time);
        fprintf('  Total time: %.1f seconds\n', simple_time + complex_time);
        
        % Update timing info
        timing_info.complex_time = complex_time;
        timing_info.total_time = simple_time + complex_time;
    else
        fprintf('Simple posterior detected - keeping initial fit results\n');
        fprintf('  Total time: %.1f seconds\n', simple_time);
    end
    
    % Store settings info in the existing nestedSamplingParams structure
    % This is a property that already exists in the class
    settings = struct(...
        'Nlive_simple', Nlive_simple, ...
        'Nlive_complex', Nlive_complex, ...
        'tolerance_simple', tolerance_simple, ...
        'tolerance_complex', tolerance_complex, ...
        'Nlive_actual', ifelse(is_complex, Nlive_complex, Nlive_simple), ...
        'tolerance_actual', ifelse(is_complex, tolerance_complex, tolerance_simple));
    
    % Merge settings and timing info
    merged_info = appendStructs(settings, timing_info);
    
    % Store in the nestedSamplingParams field which already exists in the class
    fitResult.nestedSamplingParams.twoStageInfo = merged_info;
    
    fprintf('=== Fitting Complete ===\n\n');
end

% Helper function for ternary operation
function result = ifelse(condition, true_value, false_value)
    if condition
        result = true_value;
    else
        result = false_value;
    end
end

% Helper function to merge structs
function result = appendStructs(struct1, struct2)
    result = struct1;
    fields = fieldnames(struct2);
    for i = 1:length(fields)
        result.(fields{i}) = struct2.(fields{i});
    end
end

function is_complex = analyzeComplexity(fitResult_stage1, psfInit, parEst)
    % Complexity analysis for two-stage PSF fitting
    fprintf('  Analyzing complexity...\n');
    
    complexity_count = 0;
    
    % Check 1: Parameter uncertainty analysis
    if ~isempty(fitResult_stage1.parameterUncertainties) && ~isempty(fitResult_stage1.estimatesPositionDefocus)
        uncertainties = fitResult_stage1.parameterUncertainties;
        estimates = fitResult_stage1.estimatesPositionDefocus.ML(1:end-1);
        relative_uncertainties = uncertainties ./ max(abs(estimates), 1e-6);
        max_relative_uncertainty = max(relative_uncertainties);
        
        if max_relative_uncertainty > 0.1
            complexity_count = complexity_count + 1;
            fprintf('    - High parameter uncertainties detected (%.1f%%)\n', max_relative_uncertainty * 100);
        else
            fprintf('    - Parameter uncertainties normal (%.1f%%)\n', max_relative_uncertainty * 100);
        end
    end
    
    % Check 4: Log evidence quality
    if ~isempty(fitResult_stage1.logEvidence)
        if fitResult_stage1.logEvidence < -1e99
            complexity_count = complexity_count + 1;
            fprintf('    - Very poor log evidence detected (%.1f)\n', fitResult_stage1.logEvidence);
        else
            fprintf('    - Log evidence reasonable (%.1f)\n', fitResult_stage1.logEvidence);
        end
    end
    
    % % Decision (now using 2 checks instead of 4)
    % is_complex = complexity_count >= 1;  % Changed from >= 2 to >= 1
    % Decide only based on log likelihood value
    is_complex = ~isempty(fitResult_stage1.logEvidence) && fitResult_stage1.logEvidence < -1e99;

    fprintf('    - Complexity flags: %d/2 triggered\n', complexity_count);
    if is_complex
        fprintf('    - VERDICT: Complex posterior (multimodal likely)\n');
    else
        fprintf('    - VERDICT: Simple posterior (likely Gaussian)\n');
    end
end