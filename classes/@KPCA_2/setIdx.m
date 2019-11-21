function setIdx(obj, idx)

% Idiot checks
% if length(idx) ~= size(obj.training_data, 1)
%     warning('Provided vector of cluster assignment has wrong size. Ignored.');
%     return
% end

% Set 
obj.local_idx = idx;

end


