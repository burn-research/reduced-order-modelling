function runDirectKriging(obj, varargin)


% Kriging on the training data directly
fprintf('Kriging on the original variables will be performed.\n');
obj.update_Kriging(obj.training_data, 'direct');
fprintf('Kriging terminated.\n');

% Get direct Kriging errors
obj.getDirectKrigingErrors();

end

