function [mirrored_mesh, mirrored_data, more_mirrored_data] = mirror_that(mesh, data, varargin)

if ~isempty(varargin)
    more_data = varargin{1};
else
    more_data = [];
end

mirrored_mesh = [-flip(mesh(:,1)), flip(mesh(:,2)); mesh];

if ~isempty(data)
    I_times = size(data,1) / length(mesh);
    mirrored_data = zeros(2*size(data,1), size(data,2));

    l = length(mirrored_mesh);
    for i = 1 : I_times
        dummy = data(1+length(mesh)*(i-1):i*length(mesh),:);
        mirrored_data(1+l*(i-1):i*l,:) = [flip(dummy,1); dummy];
    end
else 
    mirrored_data = [];
end

more_mirrored_data = [];
if ~isempty(more_data)
    [~, more_mirrored_data] = mirror_that(mesh, more_data);
end

end


