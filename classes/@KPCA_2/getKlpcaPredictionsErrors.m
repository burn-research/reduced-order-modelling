function getKlpcaPredictionsErrors(obj, varargin)

% If there is no test data, return
if isempty(obj.original_data)
    return
end

[obj.klpca_prediction_error_observations, ...
    obj.klpca_prediction_error_variables] = obj.getError(...
    obj.original_data, obj.klpca_predictions);

end