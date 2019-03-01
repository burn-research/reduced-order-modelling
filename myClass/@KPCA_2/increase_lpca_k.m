function obj = increase_lpca_k(obj, varargin)

k_new = varargin{1};

for i = 1 : length(obj.lpca)
    obj.lpca{i}.k = k_new;
    obj.lpca{i} = obj.lpca{i}.update();
end


end



