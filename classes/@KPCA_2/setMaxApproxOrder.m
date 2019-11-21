function setMaxApproxOrder(obj, varargin)

% Get info on size
[rows, cols] = size(obj.training_data);

% Very unlikely to find the maximum
obj.pca_approximation_order = max(min([rows, cols]) - 1, 1);

end

