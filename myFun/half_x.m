function [mesh, z] = half_x(mesh, varargin)

if ~isempty(varargin)
    z = varargin{1};
else
    z = [];
end

if length(varargin) > 1
    t = varargin{2};
else
    t = 2;
end

x_min = min(mesh(:,1));
x_max = max(mesh(:,1));
x_half = (x_max - x_min) / t;

if length(varargin) > 2
    x_half = varargin{3};
end

I = mesh(:,1) < x_half;
mesh = mesh(I,:);

if ~isempty(z)
    I_times = size(z,1) / length(I);
    Iz = repmat(I, 1, I_times);
    z = z(Iz,:);
end

end

