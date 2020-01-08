function [trimmed_modes] = trim_factor_weight(modes, threshold, iteratively)
% This function trims the factor weights. All weights below the threshold
% on a given eigenvector are set to 0.
%
% Two options exist for running this function, iteratively or not.
%
% If the non-iterative trimming is selected, the factor weights are trimmed all
% at once and the factor is re-normalized after trimming.
%
% If the iterative trimming is selected, a single weight on a given factor
% is zeroed if it is below the threshold, and the whole factor is re-normalized.
% This option makes sure that no factor will be zeroed as a result of trimming.
%
% Input:
% ------------
% - modes
%       matrix containing the factors (modes).
%
% - threshold
%       scalar specifying the trimming threshold.
%
% - iteratively
%       boolean specifying whether trimming will be performed iteratively.
%       By default it is set to false.
%
% Output:
% ------------
% - trimmed_modes
%       trimmed, normalized factors.

%% trim_factor_weight()
% Checks:
if ~exist('iteratively', 'var') || isempty(iteratively)
    iteratively = false;
end

trimmed_modes = zeros(size(modes));

if iteratively == false
        for i = 1:1:size(modes, 2)

            % Look at one column at a time:
            column = modes(:,i);

            % Set weights below threshold to zero:
            column(abs(column)<threshold) = 0;

            % Re-normalize the eigenvectors:
            trimmed_modes(:,i) = column/norm(column);

        end
else
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

end
