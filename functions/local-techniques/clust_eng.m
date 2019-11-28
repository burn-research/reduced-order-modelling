function [k] = clust_eng(X, k, n_eigs)

[rows, columns] = size(X);

% Convergence indicators initialization
convergence = 0;
iter = 0;
iter_max = 6000;
eps_rec = 1.0;
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;
r_tol = 1.0e-08;
dep = true;

% Initial time
t0 = cputime;

while dep == true

% CENTER AND SCALE THE DATA
% Select centering and scaling criteria

pre_cent_crit = 0;
pre_scal_crit = 1;

[~, X_ave] = center(X, pre_cent_crit);
[scal_X, X_gamma] = scale(X, X, pre_scal_crit);

% 1) INITIALIZATION

% Centroids initialization

C_int = linspace(1,rows, k+2);
C = scal_X(round(C_int(2:k+1)), :);

% Initialization of eigenvectors

n_eigs_max = n_eigs;
eigvec = cell(k, 1);
gamma = cell(k, 1);

for j = 1 : k
    eigvec{j} = eye(columns, n_eigs);
    gamma{j} = ones(1, columns);
end

% 2) PARTITION
% Partition can be performed using the unsupervised quantization algorithm

% Select centering and scaling criteria for PCA

cent_crit = 1;
scal_crit = 0;

while ((convergence == 0) && (iter < iter_max))
    
    C_convergence = 0;
    eps_rec_convergence = 0;   
    sq_rec_err = zeros(rows, k);
    
% Evaluate the squared mean reconstruction error
    
    for jj = 1 : k
        D = diag(gamma{jj});
        C_mat = repmat(C(jj, :), rows, 1);
        rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{jj} * eigvec{jj}' * D);
        sq_rec_err(:, jj) = sum(rec_err_os.^2, 2);
    end
        
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    
    % Evaluate the global mean error
    
    eps_rec_new = mean(rec_err_min_rel);
    [nz_X_k, nz_idx_clust, k] = partition(scal_X, idx, k);
    rec_err_min_rel_k = cell(k, 1);
    
    for jj = 1 : k
        rec_err_min_rel_k{jj} = rec_err_min_rel(nz_idx_clust{jj}, 1);
    end
    
    % Evaluate the mean error in each cluster
    
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);
    for j = 1 : k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_X_k{j}, 1);
    end
    
    % 3) EVALUATE NEW CLUSTERS' CENTROIDS
    
        C_new = zeros(k, columns);        
        for j = 1 : k
            C_new(j, :) = mean(nz_X_k{j}, 1);
        end
        eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
        fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);
        if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) ...
                && (n_eigs < n_eigs_max)) 
            n_eigs = n_eigs + 1;
            fprintf('\n Cluster %d dimension increased to %d \n', j,  n_eigs);
        end
        
% Judge convergence: clusters centroid and relative reconstruction
        % error
    
        if (eps_rec_var < r_tol)
            eps_rec_convergence = 1;
        end
        if (size(C) == size(C_new))
            C_var = abs((C_new - C) ./ (C_new + a_tol));
            if (C_var(:, :) < r_tol)
                C_convergence = 1;
            end
        end
        if ((iter > 1) && (C_convergence == 1) && (eps_rec_convergence == 1))
            convergence = 1;
            fprintf('\nConvergence reached in %d iteartion \n', iter);
        end

        % Update recontruction error and cluster centroids
        C = C_new;
        eps_rec = eps_rec_new;
        
        % 4) PERFORM LOCAL PCA
   
        % Initialization of cell arrays
        
        sort_eigval = cell(k, 1);
        sort_eigvec = cell(k, 1);
        eigval = cell(k, 1);
        eigvec = cell(k, 1);
        n_eig = cell(k, 1);
        U_scores = cell(k, 1);
        W_scores = cell(k, 1);
        gamma = cell(k, 1);
        scaled_data = cell(k, 1);
        rec_data = cell(k, 1);
       
        
        for j = 1 : k
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, rec_data{j}] ...
                = pca_lt(nz_X_k{j}, cent_crit, scal_crit, 4, n_eigs);
        end
        
        iter = iter + 1;
    
end

if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end

% ORIGINAL DATA AND RECOVERED DATA RECONSTRUCTION

rec_scal_X = zeros(rows, columns);
rec_scal_X_hat = zeros(rows, columns);
for j = 1 : k
    rec_scal_X(nz_idx_clust{j}, :) = nz_X_k{j};
    rec_scal_X_hat(nz_idx_clust{j}, :) = rec_data{j};
end
rec_cent_X_data = unscale(rec_scal_X_hat, X_gamma);
rec_X_data = uncenter(rec_cent_X_data, X_ave);

% Check that the reconstructed original data are the same of the original
% data

if (abs(scal_X(:, :) - rec_scal_X(:, :)) > a_tol)
    error('The reconstructed data are non equal to the original data');
end

% CPU TIME
overall_cpu_time = cputime - t0;

%CLUSTER ENGINEERING
for n =1:k
    A = U_scores{n};
    x = A(:,1);
    y = A(:,2);
    [pop4, gof4] = fit(x,y, 'poly4');
    if  gof4.rsquare >= 0.20 
        k = k+1;
        dep = true;
        convergence = 0;
        iter = 0;
        break
    elseif gof4.rsquare >= 0.5
        k = k+10;
        dep = true;
        convergence = 0;
        iter = 0;
        break
    elseif gof4.rsquare >= 0.3
        k = k+5;
        dep = true;
        convergence = 0;
        iter = 0;
        break
    else
        dep = false;
    end
    
    
end

fprintf('\nIterating with: %d \n', k);



end

fprintf('\nOptimized number of clusters equal to %d \n', k);
fprintf('\nIn every cluster the Scores do not exhibit any correlation %d \n');


    
% CREATE A DIRECTORY FOR THE OUTPUTS
thetime = clock;
VQ_name = 'VQPCA';
thetimestr = datestr(thetime);
dirname = [VQ_name '_n_clust_' num2str(k) '_n_eig_' num2str(n_eigs) '_' ...
    thetimestr(1:11) '_' num2str(thetime(4)) '_' num2str(thetime(5)) '_' ...
    num2str(floor(thetime(6)))];
mkdir(dirname)
cd(dirname)
fid = fopen('output.out','w');
fprintf(fid, 'Output file \n');
%fprintf(fid, '\nClusters centroids initialization: %s', opt_3);
fprintf(fid, '\n');
fprintf(fid, '\nInitial clusters centroids \n');

% WRITE THE OUTPUT FILE

fprintf(fid, '\nTotal number of clusters equal to %d \n', k);
fprintf(fid, '\nTotal number of iterations equal to %d \n', iter);
fprintf(fid, '\nRelative recontruction error equal to %d \n', eps_rec_new);
fprintf(fid, '\nRelative recontruction errors in each cluster \n', eps_rec_new_clust);
for j = 1 : k
    fprintf(fid, '%d \t', eps_rec_new_clust(j)); 
end
fprintf(fid, '\n');
fprintf(fid, '\nNumber of eigenvalues in each cluster \n');
for j = 1 : k
    fprintf(fid, '%d \t', n_eigs); 
end
fprintf(fid, '\n');
fprintf(fid, '\nTotal CPU time equal to %d s \n', overall_cpu_time);
for j = 1 : k
    fprintf(fid, '\nCluster n. %d size = ', j);
%     fprintf(fid, '%d \n', size_clust(j));
end
fprintf('\nFinal clusters centroids \n');
for j = 1 : k 
    for l = 1 : columns
    fprintf(fid, '%d \t', C(j, l));
    end
    fprintf(fid, '\n');
end

% SAVE DATA

save nz_X_k nz_X_k
save rec_X_data rec_X_data
save idx.out idx -ASCII -DOUBLE
save idx.txt idx -ASCII -DOUBLE
save C.out C -ASCII -DOUBLE
save eigval eigval
save eigvec eigvec
save U_scores U_scores

% PLOTS

    % PLOT THE PCs SCORES IN EACH LOCAL REGION
    
    for j = 1 : k
        figure,
        plot(U_scores{j}(:, 1), U_scores{j}(:,2), 'b+');
        xlabel('1st U-score');
        ylabel('2nd U-score');
        %title('PCs SCORES IN THE CLUSTER NUMBER ',{j}, ' ') 
        filetype = '.jpg';
        clust_number = num2str(j);
    end      

    %PLOT THE PCs SCORES COLORED BY CLUSTER
    
    load('U_scores.mat');
    U = cell2mat(U_scores);
    figure, scatter(U(:,1), U(:,2), 15, idx)
    
    
    %PLOT THE DISTRIBUTION OF THE ELEMENTS IN THE CLUSTER
    
    clusterpop = zeros(k, 1);
    a = numel(idx);

    for i=1:a
    for j=1:k
        if idx(i)== j
            clusterpop(j) = clusterpop(j)+1;
        end
    end
    end
    
    figure,
    bar(1:length(clusterpop), clusterpop); grid on;
    title('Cluster populations distribution')
    
end


%%%%%%%%%%%%%%%%%%%%