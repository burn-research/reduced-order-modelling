function recoverLcpca(obj, varargin)
% Initialize matrix
obj.local_recovered_data_cpca = zeros( size(obj.training_data) );
% Recover data from Local Cpca
for ii = 1 : length(obj.local_pca)
    % Recover data in one cluster
    obj.local_pca{ii}.recoverCpca(); 
    % Fill the matrix
    if ~obj.clustering_dimension
        % Write the recovered info in the correct rows
        % Write the recovered info in the correct rows
        obj.local_recovered_data_cpca(obj.local_idx == ii,:) = ...
            obj.local_pca{ii}.cpca_recovered_data;                   
    else
        % Write the recovered info in the correct columns
        obj.local_recovered_data_cpca(:,obj.local_idx == ii) = ...
            obj.local_pca{ii}.cpca_recovered_data;
    end
end
% Data need to be unscaled-uncentered
obj.local_recovered_data_cpca = obj.uncenterUnscale(... 
    obj.local_recovered_data_cpca, ...
    obj.mean_column, ...
    obj.scaling_factors);

end
