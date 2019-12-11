function [uncentered_data] = uncenter(centered_data, centerings)
% This function uncenters the centered data set.
%
% Input:
% ------------
% - centered_data
%         a centered data set.
%
% - centerings
%         a vector of centerings that were applied to the uncentered data set.
%
% Output:
% ------------
% - uncentered_data
%         an uncentered data set.

%% uncenter()
% Checks:
[~, n_vars] = size(centered_data);

if length(centerings) ~= n_vars
    error('Centerings vector must have the same number of entries as the number of variables in a data set.');
end

% Uncenter the data set:
uncentered_data = zeros(size(centered_data));

for ii = 1:1:n_vars
    uncentered_data(:,ii) = centered_data(:,ii) + centerings(ii);
end

end
