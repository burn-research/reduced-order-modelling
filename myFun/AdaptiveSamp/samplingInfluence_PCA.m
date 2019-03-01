function [Infl_rel, Infl] = samplingInfluence_PCA(X, Y, cent_crit, scal_crit, n_modes, varargin)
%% Description
%{
part of: "Adaptive Sampling for basis improvements"

Evaluation of the Influence of each snapshot: 
- on the PCA basis;
- on the regression model;

INPUTS
    X (obs x var): input locations
    Y (obs x var): data matrix
    cent_crit: PCA centering criterion (usually no centering)
    scal_crit: PCA scaling criterion (usually no scaling)

OUTPUTS
    Infl_rel: relative influence of each observation
%}

%% Input check
% The centering and scaling criteria are usually zero, although this
% function still lets the user choose.
% Centergin criterion
if ~exist('cent_crit', 'var')|| isempty(cent_crit)
    cent_crit = 0;
end
% Scaling criterion
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end
% Number of modes
if ~exist('n_modes', 'var') || isempty(n_modes)
    n_modes = min(size(Y));
end
% Verbose
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
% Number of points
n_points = size(X, 1);

%% PCA 
fprintf('[samplingInfluence_PCA] PCA... \n');
pca_model = pcaFun(Y, cent_crit, scal_crit); % Creates a PCA structure

%% Snapshot influences on the modes
a_tol = 1e-16;
% Memory allocation
pod_nosnap = cell(n_points, 1);
Infl = [];
% Create the matrix of snapshot influences on the modes
for j = 1 : n_points
    fprintf('[samplingInfluence_PCA] %i out of %i \n', j, n_points);
    % Create the matrix with the missing snapshot
    Yj = Y;   
    Yj(j,:) = []; 
    % PCA with a missing snapshot
    pod_nosnap{j} = pcaFun(Yj, cent_crit, scal_crit); 
    % Infl(i,j): influence of snapshot j on the pca mode i
    k = min([pod_nosnap{j}.k, n_modes]);
    for i = 1 : k
        modeprod = abs( pca_model.modes(:,i)' * pod_nosnap{j}.modes(:,i) );
        Infl(i, j) = 1 / (modeprod + a_tol) - 1; 
    end
end

%% Influence of the snapshots on the modal basis
% Memory allocation
Infl_s = zeros(length(pod_nosnap), 1);
% Influence of the snapshot j on the modal basis
for j = 1 : length(pod_nosnap)
    k = min([pod_nosnap{j}.k, n_modes]);
    sv = sqrt(pod_nosnap{j}.eigenv(1:k)); sv = sv(:); % Singular values
    Infl_s(j) = sv' * Infl(:,j);
end
% Relative influence of the jth snapshot on the modal basis
Infl_rel = Infl_s / sum(Infl_s);

end

% Small function for PCA
function this = pcaFun(Y, cent_crit, scal_crit)
[this.Ycs, this.m] = center(Y, cent_crit); % Y_ave is a vector
[this.Ycs, this.d] = scale(this.Ycs, Y, scal_crit); % Y_gamma is a vector
% Apply PCA (Using Matlab's default toolbox)
[this.modes, ~, this.eigenv] = pca(this.Ycs, 'algorithm', 'svd', 'centered', false);
% Approximation order
this.k = size(this.modes, 2);      
end













