function [idx, scores] = feature_assisted_clustering(X, modes, cent_crit, scal_crit)
% This function performs feature assisted clustering based on user-input
% modes.
% NOTE: DATA WILL BE TRUNCATED TO SPECIES MASS FRACTIONS ONLY. TEMPERATURE
% IS REMOVED FROM THE ANALYSIS.
%
% Outputs: ---
% 
% - idx
%       
%       clusters indices.
%
% - scores
%
%       PC-scores.
%
% Inputs: ---
%
% - X
%       data matrix, uncentered, unscaled.
%       Data will be truncated to species mass fractions only. Temperature
%       is removed from the analysis.
%
% - modes
%
%   '0' if the modes are artificially assigned,
%   otherwise it is a matrix with the modes coming from PCA, NNMF, ICA etc.
%
% - cent_crit
%       
%       centering criteria as per `center` function.
%
% - scal_crit
%
%       scaling criteria as per `scale` function.

%% feature_assisted_clustering
% Remove temperature:
    if modes == 0
X = X(:, 2:end);
    end
% Center and scale data:
[X_cent_scal, ~] = center(X, cent_crit);
[X_cent_scal, ~] = scale(X_cent_scal, X, scal_crit);

n_obs = size(X_cent_scal, 1);
idx = zeros(n_obs,1);

% Use artificial modes:
    if modes == 0
        
        group_fuel = [1,0,0,0,0,0,0,0,1,0,0,0];
        group_oxi = [0,1,0,0,0,0,0,0,0,0,0,1];
        group_light = [0,0,1,1,0,0,0,0,0,0,0,0];
        group_h = [0,0,0,0,0,1,0,0,0,0,0,0];
        group_prod = [0,0,0,0,1,0,0,0,0,1,0,0];
        group_hco = [0,0,0,0,0,0,0,0,0,0,1,0];
        group_hxox= [0,0,0,0,0,0,1,1,0,0,0,0];
        fuel = group_fuel/norm(group_fuel);
        light = group_light/norm(group_light);
        h = group_h/norm(group_h);
        prod = group_prod/norm(group_prod);
        hco = group_hco/norm(group_hco);
        hxox = group_hxox/norm(group_hxox);
        oxi = group_oxi/norm(group_oxi);

        artificial_from_features = [fuel' oxi' light' h' hxox' hco' prod'];

        scores = X_cent_scal*artificial_from_features;

        for i = 1:n_obs
            [~, I] = max(abs(scores(i,:)));
            idx(i,1) = I;
        end
        
    else

        scores = X_cent_scal*modes;
        
        for i = 1:n_obs
            [~, idx(i)] = max(scores(i,:));
        end
        
    end

end