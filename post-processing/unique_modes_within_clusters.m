function [unique_modes, no_unique] = unique_modes_within_clusters(modes_clusters, k, alpha, explain)
% This function finds how many unique modes within local clusters there are.
% The modes are compared based on a dot product. Two modes are considered
% the same if their dot product is larger than or equal to alpha.
%
% This function also prints the mode analysis on screen.
%
% NOTE: It is assumed that all modes are normalized (their dot product must
% fall between -1 and +1).
%
% INPUT PARAMETERS --------------------------------------------------------
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
% - explain
%
%     is a boolean. If set to true there will be a verbal explanation of
%     which modes repeat in which clusters.
%
% OUTPUT PARAMETERS -------------------------------------------------------
%
% - unique_modes
%
%     a matrix with unique eigenvectors.
%
% - no_unique
%
%     is the number of uniqe modes within clusters.

% Checks:
if mod(size(modes_clusters, 2), k) ~= 0
    error(['The local modes do not divide evenly between ', num2str(k), ' clusters. Change the number of clusters `k`.']);
end

if ~exist('explain')
    explain = false;
end

%% Search for unique modes

% Perform cross-check for repetitions within local modes inside clusters:
no_modes_clusters = size(modes_clusters, 2);
modes_per_cluster = no_modes_clusters / k;
no_repeating = 0;
no_unique = 0;
rep_modes_clusters = [];
cluster_annotation = [];

for i = 1:1:no_modes_clusters

    % This mode will be compared with no_modes_clusters - 1
    extracted_mode_cluster = modes_clusters(:,i);

    if explain == true

        % Get number of cluster from which the current mode has been extracted:
        in_cluster = ceil(i/(no_modes_clusters/k));

        % Get the current mode number within the cluster:
        mode_no = rem((i - ((in_cluster-1)*modes_per_cluster)), modes_per_cluster + 1);

        disp(['Mode number ', num2str(mode_no), ' from cluster ', num2str(in_cluster)]);

    end

    for j = i+1:1:no_modes_clusters

        extracted_mode_for_comparison = modes_clusters(:,j);

        if abs(dot(extracted_mode_cluster, extracted_mode_for_comparison)) >= alpha

            if explain == true

                in_cluster_found = ceil(j/(no_modes_clusters/k));
                mode_no_found = rem((j - ((in_cluster_found-1)*modes_per_cluster)), modes_per_cluster + 1);

                disp(['    is found as mode ', num2str(mode_no_found),' in cluster ', num2str(in_cluster_found)])
                disp([' '])

            end

            modes_clusters(:,j) = zeros(size(modes_clusters, 1), 1);

        end

    end

end

unique_modes = modes_clusters(:, any(modes_clusters));
no_unique = size(unique_modes, 2);

%% Verbose:
disp(['Out of ', num2str(no_modes_clusters), ' modes in local clusters, there are ', num2str(no_unique), ' unique modes'])
disp(['using alpha = ', num2str(alpha*100), '%.'])

end
