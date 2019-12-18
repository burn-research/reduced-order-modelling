function [repeating_modes_set_1, repeating_modes_set_2, no_repeating] = unique_modes(modes_set_1, modes_set_2, offset)
% This function finds the unique modes between two sets of modes. 
% The modes are compared based on a dot product. Two modes are considered
% the same if their dot product does not differ by more than the offset.
%
% INPUT PARAMETERS --------------------------------------------------------
%
% - mode_set_1, mode_set_2
%
%     a matrix with eigenvectors.
%     Each column of this matrix is an eigenvector belonging to a
%     particular set.
%
%       Eigvec-1     Eigvec-M
%     [                       ] weight 1
%     [                       ] weight 2
%     [                       ]   .
%     [       Mode set 1      ]   .
%     [                       ]   .
%     [                       ]   .
%     [                       ] weight M
%
%     It is assumed that all mode sets have the same number of weights.
%
% - offset
%
%     user-specified offset which defines a numerical value up to which
%     modes are allowed to differ.
%
% OUTPUT PARAMETERS -------------------------------------------------------
%
% - repeating_modes_set_1, repeating_modes_set_2
%
%     are the paired matrices containing modes that repeat between both
%     sets. Both are saved because they might not be exactly the same
%     vectors due to the imposed offset.
%
% - no_repeating
%
%     is the number of modes that repeat between set 1 and set 2.

% Checks:
if size(modes_set_1, 1) ~= size(modes_set_2, 1)
    error(['Different number of variables (weights) between both sets of modes.']);
end

%% Search for unique modes:
no_modes_set_1 = size(modes_set_1, 2);
no_modes_set_2 = size(modes_set_2, 2);
no_repeating = 0;
repeating_modes_set_1 = [];
repeating_modes_set_2 = [];

for i = 1:1:no_modes_set_1

    extracted_mode_set_1 = modes_set_1(:,i);

    for j = 1:1:no_modes_set_2

        extracted_mode_set_2 = modes_set_2(:,j);

        if abs(dot(extracted_mode_set_1, extracted_mode_set_2)) >= (1-offset)
            repeating_modes_set_1 = [repeating_modes_set_1, extracted_mode_set_1];
            repeating_modes_set_2 = [repeating_modes_set_2, extracted_mode_set_2];
        end

    end

end

no_repeating = size(repeating_modes_set_1, 2);
    
end