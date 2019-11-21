function [varargout] = applyConstraint(obj, varargin)
%% Description
% Apply the constraint on the data directly (in a kind of hardcoded way).
% The routine just edits the final predictions and sets them all positive,
% as this is what is usually needed, but the routine needs to be more
% general.
%

%% Default values
on_training = false; % First input
on_preditiction = true; % Second input

%% Inputs
n_args = length(varargin);
% First input
if n_args > 0
    on_training = varargin{1};
end
% Second input
if n_args > 1
    on_preditiction = varargin{2};
end

%% Sets all the recovered observations positive
if on_training
    % PCA predictions
    if ~isempty(obj.pca_recovered_data)
        I2 = obj.pca_recovered_data < 0;
        obj.pca_recovered_data(I2) = 0;
    end

    % LocalPCA+Kriging predictions
    if ~isempty(obj.local_recovered_data_pca)
        I2 = obj.local_recovered_data_pca < 0;
        obj.local_recovered_data_pca(I2) = 0;
    end
end

%% Sets all the predictions positive
if on_preditiction
    % Direct Kriging preditions
    if ~isempty(obj.kriged_direct_data)
        I2 = obj.kriged_direct_data < 0;
        obj.kriged_direct_data(I2) = 0;
    end

    % PCA+Kriging predictions
    if ~isempty(obj.kpca_predictions)
        I2 = obj.kpca_predictions < 0;
        obj.kpca_predictions(I2) = 0;
    end

    % LocalPCA+Kriging predictions
    if ~isempty(obj.klpca_predictions)
        I2 = obj.klpca_predictions < 0;
        obj.klpca_predictions(I2) = 0;
    end
end

%% Re-estimate errors
obj.getPredictionErrors();

end

