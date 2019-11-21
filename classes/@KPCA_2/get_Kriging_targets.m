function get_Kriging_targets(obj, varargin)

if ~isempty(obj.original_data)
    temp = center_scale(obj.original_data, obj.mean_column, obj.scaling_factors);
    obj.pca_scores_kriging_targets = obj.pca_modes' * temp;
end

% if ~isempty(obj.original_data) && ~obj.clustering_dimension
%     for i = 1 : length(obj.local_pca)
%         obj.local_pca{i}.original_data = obj.original_data(obj.local_nz_idx_clust{i},:);
%         obj.local_pca{i}.get_Kriging_targets();
%     end
% elseif ~isempty(obj.original_data) && obj.clustering_dimension
%     for i = 1 : length(obj.local_pca)
%         obj.local_pca{i}.original_data = obj.original_data(:,obj.local_nz_idx_clust{i});
%         obj.local_pca{i}.get_Kriging_targets();
%     end
% end

end


