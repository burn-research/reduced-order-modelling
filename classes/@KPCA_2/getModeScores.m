function [mode, scores] = getModeScores(obj, varargin)
%% Check inputs:
nin = size(varargin,2); k=1;
if nin>0
    k=varargin{1}; if ~(k==round(k)); error('2nd input: must be an integer.'); end
end


if ~isa(obj,'cell')
    error('Invalid input.');
end

%% Main:
idx = obj{1}.idx; idxsize = length(idx);
clustNum = max(idx); counter=zeros(clustNum,1);

mode = zeros(idxsize,1); a=[]; ak=[];

for l=1:idxsize
    j=idx(l); 
    counter(j)=counter(j)+1; x=counter(j);
    mode(l,1) = obj{j}.modes(x,k);
    if counter(j)<=1
        a = [a;obj{j}.a(k,:)];
        ak = [ak;obj{j}.ak(k,:)];
    end
end
scores.a = a;
scores.ak = ak;



end