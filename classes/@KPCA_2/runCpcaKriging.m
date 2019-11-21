function runCpcaKriging(obj, varargin)

% Cannot run if no training points are present
if isempty(obj.training_points)
    warning('No training points found.');
    return
end
% Kriging
n_args = length(varargin);
if n_args > 1 && ~isempty(varargin{1})
    % User-supplied kriged scores
    obj.kriged_pca_scores = varargin{1};
else
    % Kriging on the CPCA scores
    fprintf('Kriging on the CPCA scores is in progress...\n');
    temp = 'cpca';
    test_data = []; % Test data for CPCA is too costly to evaluate
    obj.cpcaKrigingModel = {};
    k = {};
    obj.mse = [];
    mse = [];
    if obj.targetBYtarget
        for i = 1 : size(obj.cpca_scores,1)
            evalc('obj.kriged_cpca_scores(i,:) = obj.update_Kriging(obj.cpca_scores(i,:), temp);');
            k{end+1} = obj.cpcaKrigingModel;
            mse = [mse, obj.mse];
        end
        obj.cpcaKrigingModel = k;
        obj.cpcaKrigingMSE = mse;
    else
        evalc('obj.update_Kriging(obj.cpca_scores, temp, test_data);');
    end
    fprintf('Kriging on the CPCA scores has terminated.\n');
end
% Get the predictions
obj.getKcpcaPredictions();
% Terminate here if this is a LPCA object
if obj.is_local
    return
end
% Estimate the errors
obj.getKcpcaPredictionsErrors();

end

