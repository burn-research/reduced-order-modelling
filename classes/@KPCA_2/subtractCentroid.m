function c_data = subtractCentroid(obj, varargin)

obj.centroid = obj.getCentroid();
if obj.is_mesh_variable
    r = obj.centroid(:)'; % To be sure it's a row-vector
    s = size(obj.training_data, 1);
    c_data = obj.training_data - repmat(r, s, 1);
else
    r = obj.centroid(:); % To be sure it's a column-vector
    s = size(obj.training_data, 2);
    c_data = obj.training_data - repmat(r, 1, s);
end

end


