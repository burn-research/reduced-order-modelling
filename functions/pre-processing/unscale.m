function [unscaled_data] = unscale(scaled_data, scalings)
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

%% unscale()
% Checks:
[~, n_vars] = size(scaled_data);

if length(scalings) ~= n_vars
    error('Scalings vector must have the same number of entries as the number of variables in a data set.');
end

% Unscale the data set:
unscaled_data = zeros(size(scaled_data));

for ii = 1:1:n_vars
    unscaled_data(:,ii) = scaled_data(:,ii) * scalings(ii);
end

end
