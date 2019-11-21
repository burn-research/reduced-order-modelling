function getKpcaPredictionsErrors_variables(obj, varargin)

obj.kpca_prediction_error_observations = ...
    obj.get_error(obj.original_data', obj.kpca_predictions');

end