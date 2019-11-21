function trainDirectKriging(obj, varargin)

% Cannot run if no training points are present
if isempty(obj.training_points)
    warning('No training points found.');
    return
end

% Train the Kriging model for Global PCA
obj.directKrigingModel = obj.train(obj.centered_scaled_data);

end
