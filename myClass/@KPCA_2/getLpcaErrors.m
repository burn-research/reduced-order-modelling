function getLpcaErrors(obj, varargin)

[obj.local_pca_reconstruction_error_observations, ...
    obj.local_pca_reconstruction_error_variables] = obj.getError(...
    obj.training_data, obj.local_recovered_data_pca);

end