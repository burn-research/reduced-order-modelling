function getKlcpcaPredictionsErrors(obj, varargin)

if isempty(obj.original_data)
    return
end

[obj.klcpca_prediction_error_observations, ...
    obj.klcpca_prediction_error_variables] = obj.getError(...
    obj.original_data, obj.klcpca_predictions);

end