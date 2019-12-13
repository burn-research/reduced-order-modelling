function [trimmed_modes] = trim_factor_weight(modes, threshold)
% This function trims the factor weights. All weights below the threshold on a given eigenvector are set to 0.
% The factor is re-normalized after trimming.
%
% Input:
% ------------
% - modes
%       matrix containing the factors (modes).
%
% - threshold
%       scalar specifying the trimming threshold.
%
% Output:
% ------------
% - trimmed_modes
%       trimmed, normalized factors.

%% trim_factor_weight()
trimmed_modes = zeros(size(modes));

    for i = 1:1:size(modes, 2)

        % Look at one column at a time:
        column = modes(:,i);

        % Set weights below threshold to zero:
        column(abs(column)<threshold) = 0;

        % Re-normalize the eigenvectors:
        trimmed_modes(:,i) = column/norm(column);

    end

end
