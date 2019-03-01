function [maxloc, pot, relInf] = AdaptiveSampling(X, Y, X_cdt, reg_based, k, a, trendFun, corrFun, verbose, varargin)
%% Description
% Adaptive Sampling for basis improvements.
% INPUTS
%     X (obs x var): starting points
%     Y (obs x var): starting values
%     X_cdt (obs x var): candidate points
%     reg_based (optional): bool, regression-based influence
%     k (optional): int, app. order for pca-based influence
% OUTPUTS
%     maxloc (1 x var): position in X_cdt of the new point to sample
%     pot: potential of each candidate point
%     relInf: relative influence of each starting point
%

%% Input check
% Number of available observations
n_obs = size(X, 1);
% Candidate points
if ~exist('X_cdt', 'var') || isempty(X_cdt)
    X_cdt = [];
    n_cdt = 0;
else
    % Be sure that no value in X is in X_cdt
    I = ismember(X_cdt, X, 'rows');
    X_cdt(I) = [];
    n_cdt = size(X_cdt, 1);
end
% Regression-based influence?
if ~exist('reg_based', 'var') || isempty(reg_based)
    reg_based = false;
end
% App. order in case of PCA-based influence
if ~exist('k', 'var') || isempty(k)
    k = min(size(Y));
end
% Correction for the distance
if ~exist('a', 'var') || isempty(a)
    a = 1.0;
end
% Trend function
if ~exist('trendFun', 'var') || isempty(trendFun)
    trendFun = 'regpoly0';
end
% Kernel
if ~exist('corrFun', 'var') || isempty(corrFun)
    corrFun = @corrmatern32;
end
% Verbose
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

%% Relative Influence
if reg_based
    % Kriging-based
    is_slow = false;
    relInf = samplingInfluence_Kriging(X, Y, trendFun, corrFun, [], is_slow);
else
    % PCA-based
    relInf = samplingInfluence_PCA(X, Y, 0, 0, k, verbose);
end

%% Enrichment Potential
if n_cdt > 0
    % Memory allocation
    pot = zeros(n_cdt, 1);
    % Potential of each candidate point
    for i = 1 : n_cdt
        dist = zeros(n_obs, 1);
        % Distance of all starting points from each candidate point
        for j = 1 : n_obs
            dist(j) = norm(X_cdt(i,:) - X(j,:));
        end
        % Find the index of the closest point (the starting point that is the
        % closest to this candidate point
        [~, minloc] = min(dist);
        % Choose the j for which it is max
        pot(i) = (dist(minloc)^a) * relInf(minloc);    
    end
    % Candidate point to be chosen: get the index of the candidate point that 
    % maximizes the Pot
    [~, maxloc] = max(pot);
else
    maxloc = [];
    pot = [];
end

end



