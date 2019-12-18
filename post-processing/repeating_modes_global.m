function [rep_modes_global, rep_modes_clusters, cluster_annotation, no_repeating] = repeating_modes_global(modes_global, modes_clusters, k, alpha)
% This function finds the repeating modes between global set of modes
% and modes in local clusters.
% The modes are compared based on a dot product. Two modes are considered
% the same if their dot product is larger than or equal to alpha.
%
% NOTE: It is assumed that all modes are normalized (their dot product must
% fall between -1 and +1).
%
% INPUT PARAMETERS --------------------------------------------------------
%
% - modes_global
%
%     a matrix with normalized eigenvectors.
%     Each column of this matrix is an eigenvector belonging to a
%     particular set.
%
%       Eigvec-1     Eigvec-M
%     [                       ] weight 1
%     [                       ] weight 2
%     [                       ]   .
%     [      Global modes     ]   .
%     [                       ]   .
%     [                       ]   .
%     [                       ] weight M
%
%     It is assumed that all mode sets have the same number of weights.
%
% - modes_clusters
%
%     a matrix with normalized eigenvectors combined from all clusters in
%     the following way: [modes_cluster_1, modes_cluster_2, ... ].
%     Each column of this matrix is an eigenvector belonging to a
%     particular cluster.
%
%       Eigvec-1     Eigvec-N
%     [                       ] weight 1
%     [                       ] weight 2
%     [                       ]   .
%     [   Modes in clusters   ]   .
%     [                       ]   .
%     [                       ]   .
%     [                       ] weight M
%
%     It is assumed that all mode sets have the same number of weights.
%
% - k
%
%     is the number of clusters represented by the matrix mode_local.
%
% - alpha
%
%     a numerical value up to which modes are allowed to differ.
%
% OUTPUT PARAMETERS -------------------------------------------------------
%
% - rep_modes_global, rep_modes_clusters
%
%     are the paired matrices containing modes that repeat between both
%     sets. Both are saved because they might not be exactly the same
%     vectors due to the imposed offset.
%
% - cluster_annotation
%
%     is a vector whose each entry is associated to a respective column of
%     rep_modes_local and specifies in which cluster this repeating
%     mode was present.
%
% - no_repeating
%
%     is the number of unique global modes that repeat in cluster modes.

% Checks:
if size(modes_global, 1) ~= size(modes_clusters, 1)
    error(['Different number of variables between both sets of modes.']);
end

if mod(size(modes_clusters, 2), k) ~= 0
    error(['The local modes do not divide evenly between ', num2str(k), ' clusters. Change the number of clusters `k`.']);
end

%% Search for unique modes

% Remove possible repetitioons from global modes matrix:
% (in case of PCA this should not happen, since all eigenvectors are
% orthogonal)
for i = 1:1:size(modes_global, 2)

    extracted_mode_global = modes_global(:,i);

    compare_global = modes_global;
    compare_global(:,i) = [];

    for j = 1:1:size(compare_global, 2)

        compare_mode_global = compare_global(:,j);

        if abs(dot(extracted_mode_global, compare_mode_global)) >= alpha
            modes_global(:,i) = zeros(size(modes_global, 1), 1);
        end

    end

end

modes_global = modes_global(:, any(modes_global));

% Perform cross-check for repetitions between global and local modes:
no_modes_global = size(modes_global, 2);
no_modes_clusters = size(modes_clusters, 2);
no_repeating = 0;
rep_modes_global = [];
rep_modes_clusters = [];
cluster_annotation = [];

for i = 1:1:no_modes_global

    extracted_mode_global = modes_global(:,i);

    for j = 1:1:no_modes_clusters

        extracted_mode_local = modes_clusters(:,j);

        if abs(dot(extracted_mode_global, extracted_mode_local)) >= alpha

            if size(rep_modes_global, 1) ~= 0
                [~, index] = ismember(extracted_mode_global', rep_modes_global', 'rows');
            else
                index = 0;
            end

            if index == 0
                rep_modes_global = [rep_modes_global, extracted_mode_global];
            end

            rep_modes_clusters = [rep_modes_clusters, extracted_mode_local];
            cluster_annotation = [cluster_annotation, ceil(j/(no_modes_clusters/k))];

        end

    end

end

no_repeating = size(rep_modes_global, 2);

%% Verbose:
disp(['Out of ', num2str(size(modes_global, 2)), ' global modes, ', num2str(no_repeating), ' are repeating in local clusters'])
disp(['using alpha = ', num2str(alpha*100), '%.'])

end
