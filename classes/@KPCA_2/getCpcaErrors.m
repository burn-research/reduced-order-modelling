function getCpcaErrors(obj, varargin)

[obj.cpca_reconstruction_error_observations, ...
    obj.cpca_reconstruction_error_variables] = obj.getError(...
    obj.training_data, obj.cpca_recovered_data);

end