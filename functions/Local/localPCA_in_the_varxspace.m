function [idx, uncentered_nz_X_k, F_k, varargout] = localPCA_in_the_varxspace(X, n_eigs)

F=X(:,1);

[rows, columns] = size(X);

% INPUTS
% Choose the number od clusters
k = input('\nSpecify the number of clusters \n');

% Convergence indicators initialization
convergence = 0;
iter = 0;
iter_max = 600;
eps_rec = 1.0;
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;
r_tol = 1.0e-08;

% CENTER AND SCALE THE DATA
% Select centering and scaling criteria
pre_cent_crit = 0;  pre_scal_crit = 1;

[cent_X, X_ave] = center(X, pre_cent_crit);
[scal_X, X_gamma] = scale(X, X, pre_scal_crit);

% 1) INITIALIZATION
% Centroids initialization
C_int = linspace(1, rows, k+2);
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
% Partition can be performed in two ways:
% 1) Iteratively, using the unsupervised quantization algorithm
% 2) By conditioning the data using a supervised quantization algorithm,
% i.e. by conditioning the data on the mixture fraction

% Select centering and scaling criteria for PCA
cent_crit = 1;  scal_crit = 0;

while ((convergence == 0) && (iter < iter_max))
    
    C_convergence = 0;
    eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);  

    sq_rec_err = zeros(rows, k);
    

    % Evaluate the squared mean reconstruction error
    for j = 1 : k
        D = sparse(diag(gamma{j}));
        C_mat = repmat(C(j, :), rows, 1);
    
    t=size(eigvec{j},2); temp=eigvec{j}; temp2=temp; %rebuilds the centered-scaled data
    d=gamma{j}; d=d(:);
    for i=1:t
        temp(:,i)=temp(:,i).*d; 
        temp2(:,i)=temp2(:,i)./(d+eps);
    end
    PI=temp2 * temp';
        rec_err_os = (scal_X - C_mat) - (scal_X - C_mat)*(PI);
        
%         rec_err_os = (scal_X - C_mat) - (scal_X - C_mat)*( D \ eigvec{j} * eigvec{j}' * D);
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2); %this is like a norm
    end
        
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    
    % Partition the data into clusters
    [nz_X_k, nz_idx_clust, k] = partitionVQ(scal_X, idx, k);

   
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
% %                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, rec_data{j}] ...
% %                = ...
% %                pca_lt_new(nz_X_k{j}, cent_crit, scal_crit, 4, n_eigs);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigval{j}, eigvec{j}, n_eig{j}, U_scores{j}, gamma{j}, scaled_data{j}, rec_data{j}
            
            % OPTIONS variables & My Class
            thisx=nz_X_k{j}';
            pcaOpts.mycon=[]; pcaOpts.con=false; pcaOpts.method=2;
            if (pcaOpts.method==2)&&(n_eigs>size(thisx,2))
                pcaOpts.k=size(thisx,2);
            elseif (pcaOpts.method==1)&&(n_eigs>size(thisx,1))
                pcaOpts.k=size(thisx,1);
            else
               pcaOpts.k=n_eigs; 
            end
            pod=pcaData(thisx,pcaOpts); % Build the pcaData object
            
            % AP's toolbox's needed variables
            eigval{j}=pod.eigenv;       eigvec{j}=pod.modes;
            n_eig{j}=pod.k;             U_scores{j}=pod.a';
            scaled_data{j}=pod.Ycs';    rec_data{j}=pcaData.Rebuild(pod)';
            gamma{j}=pod.d;
    end

    iter = iter + 1;
    
end

if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end

% MIXTURE FRACTION PARTITION
F_k = cell(k, 1);
unscaled_nz_X_k = cell(k, 1);
uncentered_nz_X_k = cell(k, 1);
for j = 1 : k
    F_k{j} = F(nz_idx_clust{j}, 1);
    unscaled_nz_X_k{j} = unscale(nz_X_k{j}, X_gamma);
    uncentered_nz_X_k{j} = uncenter(unscaled_nz_X_k{j}, X_ave(nz_idx_clust{j}, :));    
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

% VARARGOUT (ADDED BY ME)
    if (nargout>3)
        varargout{1}.nz_X_k=nz_X_k;
        varargout{1}.k=k;
%         varargout{1}.eigvec=eigvec;
%         varargout{1}.
    end

close all hidden
% cd ..


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
                
    
    