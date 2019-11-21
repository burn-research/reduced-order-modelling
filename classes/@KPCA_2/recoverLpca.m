function recoverLpca(obj, varargin)

% Initialize matrix
obj.local_recovered_centered_data = zeros( size(obj.training_data) );

% Recover data from Local Pca
for ii = 1 : length(obj.local_pca)
    % Recover data in one cluster
    obj.local_pca{ii}.recoverPca(); 
    % Fill the matrix
    if ~obj.clustering_dimension
        % Write the recovered info in the correct rows
        obj.local_recovered_centered_data(obj.local_idx == ii, :) = ...
                                obj.local_pca{ii}.pca_recovered_data;                   
    else
        % Write the recovered info in the correct columns
        obj.local_recovered_centered_data(:, obj.local_idx == ii) = ...
                                obj.local_pca{ii}.pca_recovered_data;
    end
end

% This data is still centered-scaled. We have recovered the centered-scaled
% data that was passed as input to the clustering/local PCA procedure.

obj.local_recovered_data_pca = obj.uncenterUnscale(... 
    obj.local_recovered_centered_data, ...
    obj.mean_column, ...
    obj.scaling_factors);

end



