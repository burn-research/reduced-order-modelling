function getKcpcaPredictionsErrors(obj, varargin)

if isempty(obj.original_data)
    return
end

[obj.kcpca_prediction_error_observations, ...
    obj.kcpca_prediction_error_variables] = obj.getError(...
    obj.original_data, obj.kcpca_predictions);

end