function [idx, clusters, u_scores, eigvec, nz_idx_clust, info] = lpca_partition(X, n_eigs, k, cent_crit, scal_crit, r2_min, varargin)
%% Description:
% The objective of this function is to return clusters that have low PCA
% reconstruction error and also no structure.
%

%% Input
if nargin < 4 || isempty(cent_crit)
    cent_crit = 1;
end
if nargin < 5 || isempty(scal_crit)
    scal_crit = 0;
end
if nargin < 6
    r2_min = .12;
end

%% Main
iter = 0;
iter_max = 200;
cond = true;
dep = false;
fprintf('Running Local PCA clustering and checking no structure is present. \n');
while cond
    fprintf('Current number of clusters: %i \n', k);
    evalc('[idx, clusters, u_scores, ~, ~] = localPCA(X, n_eigs, k, 1, false, cent_crit, scal_crit);');
    % Cluster engineering
    k_new = length(u_scores);
    k = k_new;
    for ii = 1 : k_new
        x = u_scores{ii}(:,1);
        y = u_scores{ii}(:,2);
        [~, gof_1] = fit(x, y, 'poly4');
        [~, gof_2] = fit(y, x, 'poly4');
        if gof_1.rsquare >= 0.5 || gof_2.rsquare >= 0.5
            k = k + 10;
            break
        elseif gof_1.rsquare >= 0.3 || gof_2.rsquare >= 0.3
            k = k + 5;
            break
        elseif gof_1.rsquare >= .2 || gof_2.rsquare >= .2
            k = k + 1;
            break
        else
            dep = true;
        end 
    end
    % Update iter
    iter = iter + 1;
    % Check convergence
    cond = ~((iter > iter_max) || dep);
end
% Re-evaluate gof4
k_new = length(u_scores);
info.gof_1 = cell(k_new, 1);
info.gof_2 = cell(k_new, 1);
for ii = 1 : k_new
    x = u_scores{ii}(:,1);
    y = u_scores{ii}(:,2);
    [~, info.gof_1{ii}] = fit(x, y, 'poly4');
    [~, info.gof_2{ii}] = fit(y, x, 'poly4');
end
fprintf('Structures eliminated with %i clusters. \n', k_new);
% Bisection section
idx_clust = cell(k_new, 1);
k_bi = 2;
fprintf('Bisection part. \n');
for ii = 1 : k_new
    % Check the target has been met by cluster ii
    if info.gof_1{ii}.rsquare >= r2_min || info.gof_2{ii}.rsquare >= r2_min
        fprintf('Splitting cluster %i \n', ii);
        evalc('[idx_clust{ii}, ~] = localPCA(clusters{ii}, n_eigs, k_bi);');
        % Get the rows of this cluster
        mask = (idx == ii);
        idx(mask) = idx(mask) + get_b(max(idx), idx_clust{ii});
    end
end
% Get clusters and run Local PCA
[clusters] = get_clusters(X, idx);
k_new = length(clusters);
fprintf('Gathering info and returning %i clusters. \n', k_new);
[eigvec, ~, ~, u_scores, ~, ~] = lpca(clusters, n_eigs, cent_crit, scal_crit);
nz_idx_clust = cell(k_new, 1);
for ii = 1 : k_new 
    nz_idx_clust{ii} = find(idx == ii);
end
% Re-evaluate gof4
info.pop_1 = cell(k_new, 1);
info.pop_2 = cell(k_new, 1);
info.gof_1 = cell(k_new, 1);
info.gof_2 = cell(k_new, 1);
try
    for ii = 1 : k_new
        x = u_scores{ii}(:,1);
        y = u_scores{ii}(:,2);
        [info.pop_1{ii}, info.gof_1{ii}] = fit(x, y, 'poly4');
        [info.pop_2{ii}, info.gof_2{ii}] = fit(y, x, 'poly4');
    end
catch
    
end

end


function b_new = get_b(n_clust_old, b)
b_new = n_clust_old + (b-1);
b_new(b == 1) = 0;
end

