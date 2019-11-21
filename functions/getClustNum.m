% Given a point 'p' (in the parameter-space), this function tells you to
% which cluster it belongs.

function clustNum = getClustNum(p,clustPCA,varargin)
% INPUTS: - p: point in the parameter-space;
%         - clustPCA: cell arrays of 'pcaData' object;
%         - facultative, true or false,

%% Input checks:
tt=1;
if ~isempty(varargin)
    tt=varargin{1};
end
if (tt~=1)&&(tt~=0)
    error('Third input must be logical.');
end

temp = [];
for j=1:length(clustPCA)
    temp=[temp;clustPCA{j}.xp];
    if (tt)
        temp=[temp;clustPCA{j}.xpk];
    end
end

if isa(temp,'double') && (length(unique(temp)) ~= length(temp)) 
    error('The fields "clustPCA{i}.xp" and "clustPCA{i}.xpk" cannot repeat the same value.');
end

%% Main:
for j=1:length(clustPCA)
    temp = clustPCA{j}.xp;
    if (tt)
        temp = [temp;clustPCA{j}.xpk];
    end
    if ~isa(p,'cell')
        if (ismember(p,temp))
            clustNum = j;
        end
    else
        if (ismember2(p,temp))
            clustNum = j;
        end        
    end
end

end

function output = ismember2(p,t)

[r,c] = size(t); [~,cp] = size(p);

if c ~= cp
    error('Number of cols must be the same!');
end

output = false;
for i=1:r
    for j=1:c
        pp = p{1,j};    tt = t{i,j};
        if (~isa(pp,'char'))&&(~isa(pp,'string'))
            l(1,j) = (pp==tt);
        else
            l(1,j) = strcmp(pp,tt);
        end
    end
    if (l==ones(1,c))
        output = true;
    end
end



end

