function [Y_lpca, Y_rec_lpca] = lpca_rec(obj)

% Check there are clusters
temp = length(obj.lpca);
if temp == 0
    fprintf('No clusters found. LPCA not applied yet probably.\n\n');
    Y_lpca = [];
    return;
end


% Recover the data
if obj.clustering_dim == 1
   % Rebuild according to obj.idx
    [Y_lpca, Y_rec_lpca] = lpca_rec_1(obj);
    
elseif obj.clustering_dim == 2
   % Rebuild according to obj.xp_sorted
    [Y_lpca, Y_rec_lpca] = lpca_rec_2(obj); 
    
else
    error('Clustering dimension is neither 1 nor 2.');
end


% Uncenter, unscale
Y_lpca = uncenterunscale(obj, Y_lpca);


end




% Uncenter, unscale
function Y_lpca = uncenterunscale(obj, Y_lpca)

for i = 1 : size(Y_lpca, 2)
    Y_lpca(:,i) = Y_lpca(:,i) .* obj.d + obj.m;
end

end



% Rebuild according to obj.idx (CASE 1)
function [A, A_pureLpca] = lpca_rec_1(obj)

% Counter for each cluster (to get the correct row)
counter = ones( 1, length(obj.lpca) );

% Return this if there is no KLPCA, only LPCA
if isempty(obj.xp_kriged)
    A = [];
end


% Initialize matrices    
A = zeros(obj.rows, obj.cols + size(obj.xp_kriged, 1));     
A_pureLpca = zeros(obj.rows, obj.cols);


% Load the KLPCA data into one matrix
for i = 1 : length(obj.lpca)
    A(obj.nz_idx_clust{i}, :) = obj.lpca{i}.Y_sorted;
    A_pureLpca(obj.nz_idx_clust{i}, :) = obj.lpca{i}.Y_rec;
end


% Initialize matrices
% A = zeros(obj.rows, obj.cols + size(obj.xp_kriged, 1));     
% A_pureLpca = zeros(obj.rows, obj.cols);
% 
% Fill the matrices
% for i = 1 : obj.rows
%     % Get the cluster that has this row
%     j = obj.idx(i);
%     
%     % Count the row of the cluster that we need to copy
%     this_row = counter(j);
%     
%     % Copy the row (KLPCA)
%     if ~isempty(obj.xp_kriged)
%         A(i,:) = obj.lpca{j}.Y_sorted( this_row, :);
%     end
%     
%     % Copy the row (pure LPCA)
%     A_pureLpca(i,:) = obj.lpca{j}.Y_rec( this_row, :);
%     
%     % Update counter in the position of the cluster just called
%     counter(j) = counter(j) + 1;
% end


end


% Rebuild according to obj.xp_sorted (CASE 2)
function [A, A_pureLpca] = lpca_rec_2(obj)

[Y_s, xp_s] = sortKriged(obj.lpca{1}.xp_sorted, obj.lpca{2}.xp_sorted, obj.lpca{1}.Y_sorted, obj.lpca{2}.Y_sorted);

num_clusts = length(obj.lpca);
if num_clusts == 2
    A = Y_s; 
    A_pureLpca = [];
    return; % You don't need to go on
end

for i = 3 : num_clusts
    [Y_s, xp_s] = sortKriged(xp_s, obj.lpca{i}.xp_sorted, Y_s, obj.lpca{i}.Y_sorted);
end

A = Y_s;
A_pureLpca = [];

end



