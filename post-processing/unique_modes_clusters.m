function [uni_modes_clusters, mode_annotation, cluster_annotation] = unique_modes_clusters(modes_global, modes_clusters, k, alpha)
% This function finds unique modes from local clusters that are not
% present within the global modes.
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
% - uni_modes_clusters
%
%     is a matrix of unique modes from local clusters that do not repeat
%     within the global modes.
%
% - mode_annotation
%
%     is a vector whose each entry is associated to a respective column of
%     uni_modes_clusters and specifies which mode number this column is.
%
% - cluster_annotation
%
%     is a vector whose each entry is associated to a respective column of
%     uni_modes_clusters and specifies in which cluster this unique
%     mode was present.

%% Checks:
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
            modes_global(:,i) = zeros(size(modes_global, 1), 1); % zeroing repetitions
        end

    end

end

modes_global = modes_global(:, any(modes_global)); % removing zeroed repetitions

% Perform cross-check for repetitions between global and local modes:
uni_modes_clusters = [];
cluster_annotation = [];
mode_annotation = [];

for i = 1:1:size(modes_clusters, 2)

    for j = 1:1:size(modes_global, 2)

        if abs(dot(modes_global(:,j), modes_clusters(:,i))) < alpha
            flag = 1; % flag is =1 whenever the local mode does not match a global one
        else
            flag = 0;
            break; % break the for loop if at least one repetition within global modes is found
        end
    end

    % If this mode was never repeated (flag is still =1), append to unique modes in clusters:
    if flag == 1
        uni_modes_clusters = [uni_modes_clusters, modes_clusters(:,i)];
        if mod(i, size(modes_clusters, 2)/k) == 0
            mode_number = size(modes_clusters, 2)/k;
        else
            mode_number = mod(i, size(modes_clusters, 2)/k);
        end
        mode_annotation = [mode_annotation, mode_number];
        cluster_annotation = [cluster_annotation, ceil(i/(size(modes_clusters, 2)/k))];
    end

end

%% Verbose:
disp(['Using alpha = ', num2str(alpha*100), '%.'])
disp(['From ', num2str(size(modes_clusters, 2)), ' local modes, ' num2str(size(uni_modes_clusters, 2)), ' are not repeating within the global modes.'])

end
