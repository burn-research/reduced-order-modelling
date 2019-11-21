function [pv, sig] = findPV(clusters, nPC, modes, eigenvals, vary, varargin)
%% Description
% Finds Principal Variables in a dataset. Applies PCA to the data matrix
% (input can be a cell array with more than one data matrix), returns list
% of PVs (ordered in descending order of importance) based on the loadings.

%% Input
% Optional arguments
if ~iscell(clusters)
    clusters = {clusters};
    modes = {modes};
    eigenvals = {eigenvals};
    vary = {vary};
end

%% Main
n_clust = length(clusters);
pv = cell(n_clust,1);
sig = cell(n_clust,1);
for ii = 1 : n_clust
    [pv{ii}, sig{ii}] = findPV_local(clusters{ii}, nPC, modes{ii}, eigenvals{ii}, vary{ii});
end
end

function [ikeep, var_tot] = findPV_local(Y, nPC, modes, eigenvals, vary)
%% Input
% Optional arguments
y_dim = size(Y,2);
cond = true;
neta = y_dim;
while(cond)
    neta = neta - 1; % # of variables to discard
    ikeep = principal_variables(Y, neta, modes);
    var_tot = sum(vary(ikeep));
    if sum(eigenvals(1:nPC)) <= var_tot
        cond = false;
    end
end

end



