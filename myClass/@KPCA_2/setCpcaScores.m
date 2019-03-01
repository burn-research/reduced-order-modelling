function setCpcaScores(obj, varargin)

nin = length(varargin);

if nin < 2 || ~varargin{2}
    obj.cpca_scores = varargin{1};
else
    obj.kriged_cpca_scores = varargin{1};
end
    
end



