function [samples, values] = memory_test(samples, values, imv, size_m1, size_m2, varargin)

smp = [];
if ~isempty(varargin)
    smp = varargin{1};
end

% In case of memory problems
n_samples = size(samples, 1); % Number of training samples
try
    % No memory problems
    temp = zeros(n_samples, n_samples);
    temp = [];
catch ME
    % Memory problems
    warning('Because of memory issues, less samples will be used.');
    % We'll find a subset of the training samples to use
    n_total = n_samples; % Initial total number of samples
    stop_rule = true; % Flag
    while(stop_rule) 
        n_total = round(n_total * 0.8); % Decrease the current number of samples
        % Find a subset of the training data
        [samples_new, values_new] = selectTrainingSet(...
            imv, size_m2, size_m1, samples, values, n_total, ...
            smp);
        % Test that there are no more memory problems: add 10% to be sure
        n_total = size(samples_new, 1);
        n_total = round(n_total * 1.10);
        try
            % Check memory will be fine
            temp = zeros(n_total, n_total);
            stop_rule = false; 
            temp = []; 
            % These will thus be the new training dataset
            samples = samples_new;
            values = values_new;
        catch
            temp = [];
            stop_rule = true; % Try again
            warning('Memory issues still present. n_total = %d', n_total);
        end
    end
    warning('A subset of grid points was chosen. n_total = %d', n_total);
end

end


