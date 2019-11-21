function getLcpcaErrors(obj, varargin)

[obj.local_cpca_reconstruction_error_observations, ...
    obj.local_cpca_reconstruction_error_variables] = obj.getError(...
    obj.training_data, obj.local_recovered_data_cpca);

end