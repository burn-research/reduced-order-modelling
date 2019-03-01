function [varargout] = getClusterPopulations(obj, varargin)
%% Description
% Get population for each cluster. This will be achieved by using the
% property 'local_idx' only.
%

%% Pre-processing
if ~isempty(obj.local_idx)
    n_clusters = numel(unique(obj.local_idx));
    if n_clusters ~= max(obj.local_idx)
        warning('Max integer in LOCAL_IDX appears not to be the actual number of clusters.');
    end
else
    n_clusters = 1;
end

if n_clusters < 2
    fprintf('\nThere`s only 1 cluster.\n');
    return
end

%% Main
pop = zeros(n_clusters, 1);
parfor i = 1 : n_clusters
    I = (obj.local_idx == i);
    pop(i) = numel(obj.local_idx(I));
end

%% Output
if nargout > 0
    varargout{1} = pop;
end

end


