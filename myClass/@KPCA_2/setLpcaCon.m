function obj = setLpcaCon(obj, varargin)

    nin = length(varargin);
    c = true;
    if nin > 0
        c = varargin{1};
    end
    if ~islogical(c)
        c = true;
    end


    for i = 1 : obj.numberOfClusters
        obj.lpca{i}.con = c;
        obj.lpca{i} = obj.lpca{i}.update();
    end


end

