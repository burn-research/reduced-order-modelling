function cleanLpca(obj, varargin)

obj.local_idx = [];
obj.local_pca = {};
obj.local_recovered_centered_data = [];
obj.local_recovered_data_pca = []; 
obj.local_pca_reconstruction_error_observations = [];
obj.local_pca_reconstruction_error_variables = []; 
obj.local_cpca_reconstruction_error_observations = [];
obj.local_cpca_reconstruction_error_variables = [];
obj.klpca_predictions = [];
obj.klpca_prediction_error_observations = [];
obj.klpca_prediction_error_variables = [];

end

