function runLcpca(obj, varargin)
% App order
k = obj.pca_approximation_order;
% Logical conditions
cond1 = (length(obj.local_idx) == size(obj.centered_scaled_data,1));
cond2 = (length(obj.local_idx) == size(obj.centered_scaled_data,2));
cond = obj.clustering_dimension;
% Perform LCPCA
for ii = 1 : length(obj.local_pca)
    obj.local_pca{ii}.my_constraint = obj.my_constraint;
    % Pass the correct data to be approximated by CPCA
    if cond1 && ~cond
        % Dimensions coherent with ~clustering_dimension
        obj.local_pca{ii}.runCpca(k, obj.centered_scaled_data(obj.local_idx == ii,:));
    elseif cond2 && cond
        % Dimensions coherent with clustering_dimension
        obj.local_pca{ii}.runCpca(k, obj.centered_scaled_data(:,obj.local_idx == ii));
    elseif cond1
        % Flags failed, check only the dimensions
        obj.local_pca{ii}.runCpca(k, obj.centered_scaled_data(obj.local_idx == ii,:));
    else
        % Flags failed, check only the dimensions
        obj.local_pca{ii}.runCpca(k, obj.centered_scaled_data(:,obj.local_idx == ii));
    end
end

% Recover data
obj.recoverLcpca();

% Estimate errors
obj.getLcpcaErrors();


end



