function obj = updateRebuiltAndSortedData(obj, varargin)



% Recovered data PCA
obj.Y_rec = KPCA.Rebuild(obj, false, false, obj.k); % Kriged = false; Let it be CPCA in case it is


% Recovered data from the A_KRIGED
obj.Y_kriged = KPCA.Rebuild(obj, true, false, obj.k);


% Sorted Kriged and PCA data
[obj.Y_sorted, obj.xp_sorted] = sortKriged(obj.xp, obj.xp_kriged, obj.Y_rec, obj.Y_kriged);


% Sorted PCA scores
[obj.a_sorted, ~] = sortKriged(obj.xp, obj.xp_kriged, obj.a, obj.a_kriged);


% K-LPCA recovered data
obj.Y_lpca = lpca_rec(obj);


% Recovered data from LPCA
obj.Y_rec_lpca = lpca_rec(obj);



end




