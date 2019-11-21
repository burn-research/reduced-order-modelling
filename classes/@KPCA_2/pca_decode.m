function [decoded_data, varargout] = pca_decode(obj, modes, scores, scaling_factors, mean_column, varargin)

n_args = length(varargin);
if n_args > 0
    decoded_data = pca_decode(modes, scores, scaling_factors, ...
        mean_column, obj.is_mesh_variable, obj.is_local, varargin{1});
else
    decoded_data = pca_decode(modes, scores, scaling_factors, ...
        mean_column, obj.is_mesh_variable, obj.is_local);
end

end


