function [trimmed_modes] = trim_factor_weight_iteratively(modes, threshold)
% This function trims the factors weights iteratively.
% Whenever the weight on a given factor is below the threshold,
% it is set to 0 and the whole factor is re-normalized.
% This function makes sure that no factor will be zeroed.
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

%% trim_factor_weight_iteratively()
trimmed_modes = zeros(size(modes));

for i = 1:1:size(modes, 2)

    % Look at one column at a time:
    column = modes(:,i);

    % This iteration assures that any vector will not be zero:
    for j = 1:1:length(column)

        % Find the lowest in the column:
        column(column == 0) = NaN;
        [lowest_abs, index_abs] = min(abs(column));
        column(isnan(column)) = 0;

        if lowest_abs < threshold

            % Set weights below threshold to zero:
            column(index_abs, 1) = 0;

            % Re-normalize the eigenvectors:
            column = column/norm(column);

        end

    end

    trimmed_modes(:,i) = column;

end

end
