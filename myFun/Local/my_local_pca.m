 % LOCAL PCA - Performs local pca on a data set.
% 
% function [idx, uncentered_nz_X_k, F_k] = local_pca(X, F, names, C)
% INPUTS
% 
% X            = Matrix of state variables
% F            = Mixture fraction
% names        = State variables names (MUST be cellstr format)
% OPTIONAL:
% C            = cluster centroids initialization
%
% OUTPUTS
%
% idx               = Vector of indexes denoting the cluster number for
%                     each of the original data points
% uncentered_nz_X_k = Cell matrix of the partitioned data, i.e. each  cell
%                     is a cluster of points
% F_k               = Cell matrix of mixture storing the values of mixture
%                     fraction in each bin identified by LPCA
%
% This routine peforms local PCA on a data set. Two different partitions
% algorithms can be selected: a supervised (FPCA) and unsupervised (VQPCA)
% one. 
% 1) When using FPCA data are partitioned into bins of mixture fraction
% and then PCA is performed locally in each one of the bin. 
% 2) When VQPCA is performed, the data are assigned into clusters based on
% their reconstruction distance, i.e. the diference between the
% actual p-dimesnional point and its q-dimesnional approximation. The
% reconstruction error is then defined as eGSRE=X-Xq. The procedure is then
% iterative and it follows the following steps: 
% (a) Initialization: the cluster centroids are randomly chosen from the
% data set and the covariance matrix is initialized to the identity matrix
% for each cluster
% b) Partition: each observation from the sample is assigned to a cluster
% using the squared reconstruction distance given by eGSRE
% c) Update: the clusters� centroids are updated on the basis of
% partitioning
% d) Local PCA: PCA is performed in each disjoint region of the sample.
% e) Steps 2-4 are iterated until convergence is reached.
% The goodness of reconstruction given by VQPCA is measured with respect to
% the mean variance in the data as eGSRE,n=E(eGSRE)/E(var(X)). If the auto
% scaling criterion is adopted the mean variance of X after scaling equals
% one and the reconstruction error metric becomes eGSRE,n=E(eGSRE).
% Convergence can be judged using the following criteria:
% a) The normalized global scaled reconstruction error, eGSRE,n, is below a
% specific threshold, eGSRE,n.
% b) The relative change in clusters� centroids between two successive
% iterations is below a fixed threshold, i.e. 10-8.
% c) The relative change in eGSRE,n between two successive iterations is
% below a fixed threshold, i.e. 10-8. Requirements b) and c) are
% particularly useful if an explanatory analysis on the performances of
% VQPCA in terms of eGSRE,n is of interest. In this case, requirement a)
% can be relaxed and the variation of eGSRE,n as a function of the number
% of eigenvalues and clusters can be analyzed by enforcing requirements b)
% and c). Otherwise, all the three conditions can be used and an iterative
% procedure for the determination of the number of eigenvalues required to
% achieve a fixed eGSRE,n could be employed. Staring with q=1, the number
% of eigenvalues can be increased progressively until the desired error
% level is reached.

function [idx, varargout] = my_local_pca(X, k, n_eigs)

[rows, columns] = size(X);

% INPUTS
% Choose the number od clusters
% k = input('\nSpecify the number of clusters \n');
plots = 0;% input('\nDo you want to plot the results? [1/0] ');
% n_eigs = input('\nChoose the number of eigenvalues n_eigs = ');

% Set font size
set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);
set(0, 'defaulttextfontweight', 'Bold')

% Convergence indicators initialization
convergence = 0;
iter = 0;
iter_max = 600;
eps_rec = 1.0;
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;
r_tol = 1.0e-08;

% Initial time
t0 = cputime;

% CENTER AND SCALE THE DATA
% Select centering and scaling criteria
pre_cent_crit = 1;
pre_scal_crit = 1;

[cent_X, X_ave] = center(X, pre_cent_crit);
[scal_X, X_gamma] = scale(X, X, pre_scal_crit);

% Choose the clustering criterion
VQ = 1;% input('\nChoose the partition algorithm: 1) VQ-PCA 2) F-PCA \n');
if VQ == 1
    VQ_name = 'VQPCA';
else
    VQ_name = 'FPCA';
end

% 1) INITIALIZATION

% Centroids initialization

if VQ == 1
    if nargin == 3
        fprintf('\nClusters centroids not specified, initialization required \n');
        init = 1;% input('\nChoose the initialization criterion: 1) Uniform 2) Random 3) From FPCA \n');
        if init == 1
            C_int = linspace(1, rows, k+2);
            C = scal_X(round(C_int(2:k+1)), :);
            opt_3 = 'uniform';
        elseif init == 2
            C_int = randomsample(rows, k);
            C = scal_X(C_int, :);
            opt_3 = 'random';
        elseif init == 3
            opt_3 = 'from FPCA';
        end
    else
        fprintf('\nClusters centroids initialized \n');
        opt_3 = 'Initialized';
    end
elseif VQ == 2
    opt_3 = 'Data partitioned on F';
end

% Initialization of eigenvectors
n_eigs_max = n_eigs;
eigvec = cell(k, 1);
gamma = cell(k, 1);

for j = 1 : k
    eigvec{j} = eye(columns, n_eigs);
    gamma{j} = ones(1, columns);
end

% CREATE A DIRECTORY FOR THE OUTPUTS
% thetime = clock;
% thetimestr = datestr(thetime);
% dirname = [VQ_name '_n_clust_' num2str(k) '_n_eig_' num2str(n_eigs) '_' ...
%     thetimestr(1:11) '_' num2str(thetime(4)) '_' num2str(thetime(5)) '_' ...
%     num2str(floor(thetime(6)))];
% mkdir(dirname)
% cd(dirname)
% fid = fopen('output.out','w');
% fprintf(fid, 'Output file \n');
% fprintf(fid, '\nClusters centroids initialization: %s', opt_3);
% fprintf(fid, '\n');
% fprintf(fid, '\nInitial clusters centroids \n');
% if (VQ == 1 && init ~= 3)
%     for j = 1 : k 
%         for l = 1 : columns
%             fprintf(fid,'%d \t', C(j, l));
%         end
%         fprintf(fid,'\n');
%     end
% end

% 2) PARTITION

% Partition can be performed in two ways:
% 1) Iteratively, using the unsupervised quantization algorithm
% 2) By conditioning the data using a supervised quantization algorithm,
% i.e. by conditioning the data on the mixture fraction

% Select centering and scaling criteria for PCA
cent_crit = 1;
scal_crit = 0;

while ((convergence == 0) && (iter < iter_max))
    
    C_convergence = 0;
    eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);  

    sq_rec_err = zeros(rows, k);
    
    if ((VQ == 2 || init == 3) && iter == 0) 
        F_min = min(F);
        F_max = max(F);
        %F_stoich = input('\nInput the value of the stoichiometric mixture fraction: Fstoich = \n');
        %F_stoich = 0.4375;% DNS CO/H2 Flame
        %F_stoich = 0.351; % Flame F
        F_stoich = 0.0579; %JHC Flame
        [nz_X_k, nz_idx_clust] = condition(scal_X, F, k, F_min, F_max, F_stoich);
        C = zeros(k, columns);
        for j = 1 : k
            C(j, :) = mean(nz_X_k{j}, 1);
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, ...
                rec_data{j}] = pca_lt(nz_X_k{j}, cent_crit, scal_crit, 4, n_eigs);
            for l = 1 : columns
                fprintf(fid,'%d \t', C(j, l));
            end
            fprintf(fid,'\n');
        end

    end

    % Evaluate the squared mean reconstruction error
    for j = 1 : k
        D = diag(gamma{j});
        C_mat = repmat(C(j, :), rows, 1);
        rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2);
    end
        
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    
    if VQ == 1
        % Partition the data into clusters
        [nz_X_k, nz_idx_clust, k] = partition(scal_X, idx, k);
    end
   
    fprintf('\nClusters dimension \n');
    disp(nz_X_k);
    
    % Evaluate the relative recontruction errors in each cluster
    rec_err_min_rel_k = cell(k, 1);
    for j = 1 : k
        rec_err_min_rel_k{j} = rec_err_min_rel(nz_idx_clust{j}, 1);
    end
    
    % Evaluate the mean error in each cluster
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);
    for j = 1 : k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_X_k{j}, 1);
    end
    fprintf('\nGlobal mean recontruction error at iteration n. %d equal to %d \n', iter, eps_rec_new);
    fprintf('\nLocal mean recontruction error at iteration n. %d \n', iter);
    for j = 1 : k
        fprintf('%d \t', eps_rec_new_clust(j));
    end
    fprintf('\n');
    
    if VQ == 1
        
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
    elseif VQ == 2
        convergence = 1;    
    end
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





% VARARGOUT (ADDED BY ME)
    if (nargout>3)
        varargout{1}.nz_X_k=nz_X_k;
        varargout{1}.k=k;
%         varargout{1}.eigvec=eigvec;
%         varargout{1}.
    end

end                
    
    