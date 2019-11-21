function getNumberOfSamples(obj, varargin)

nSamples = size(obj.p_orig, 1);

obj.a = .65;
nin = length(varargin);
if nin > 0
    obj.a = varargin{1};
end

obj.nSamplesMax = round(obj.a * nSamples);   
    
end