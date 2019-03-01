function centroid = getCentroid(obj, varargin)

centroid = mean(obj.training_data, 1) * 0;
if ~obj.is_mesh_variable
    centroid = mean(obj.training_data, 2);
end

end

