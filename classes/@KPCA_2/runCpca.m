function runCpca(obj, varargin)
%% Inputs
n_args = length(varargin);
if n_args > 0
    % Run CPCA with this approximation order
    if varargin{1} > 0
        obj.pca_approximation_order = varargin{1};
    end
end
% Data to approximate
if n_args > 1
    min_data = varargin{2};
else
    min_data = obj.centered_scaled_data;
end

%% CPCA
if ~obj.cpca_all
    % Evaluate the CPCA scores for one value of pca_approximation_order
    obj.update_gamma(min_data);
else
    % Evaluate the CPCA scores for all values of pca_approximation_order
    for ii = 1 : min(size(obj.training_data))
        obj.pca_approximation_order = ii;
        obj.update_gamma(min_data);
        obj.cpca_scores_stored{ii} = obj.cpca_scores;
    end
end

% Recover data
obj.recoverCpca();

% Estimate errors
obj.getCpcaErrors();

end





