function getKlcpcaPredictions(obj, varargin)
a_tol = 1e-16;
% Specific cluster
idx_clust = [];
if ~isempty(varargin)
    idx_clust = varargin{1};
end

% Choose according to the clustered dimension of the data matrix
if ~obj.clustering_dimension
    % Work differently for input varargin{1} 
    if ~isempty(idx_clust)
        % If a specifi cluster is indicated, only work for that
        obj.local_pca{idx_clust}.getKcpcaPredictions();
        I = obj.local_nz_idx_clust{idx_clust};
        obj.klcpca_predictions(I,:) = obj.local_pca{idx_clust}.kcpca_predictions;
    else
        % Rebuild from all clusters
        for ii = 1 : length(obj.local_pca)
            obj.local_pca{ii}.getKcpcaPredictions();
            I = (obj.local_idx == ii);
            obj.klcpca_predictions(I,:) = obj.local_pca{ii}.kcpca_predictions;
        end
    end
else
    % Get useful sizes
    n_predictions = size(obj.prediction_points, 1);
    n_clusters = length(obj.local_pca);
    % Rebuild from all clusters
    dist = zeros(n_clusters, n_predictions);
    for ii = 1 : n_clusters
        obj.local_pca{ii}.getKxpcaPredictions();
        % Get distance of each cluster from every prediction point
        d = getDistances(obj.local_pca{ii}.training_points, obj.prediction_points, obj.dist_on_cs_points);
        dist(ii,:) = d(:)';
    end
    % Rebuild for each prediction point
%     obj.klpca_predictions = combinePredictions(predictions, mse, kriged);
    dim_out = size(obj.local_pca{1}.kxpca_predictions, 1);
    for jj = 1 : n_predictions
        % Weights evaluated from the meta-features
%         weights = 1 / dist(:,jj); % Get weights for each prediction
%         tot_w = sum(weights); % Sum of the weights
        % For each cluster, get the j-th prediction and weigh it
        % accordingly. The prediction for the j-th prediction point will be
        % a weighed average of the clusters predictions, with the distance
        % evaluated before as weights.
        % For prediction jj, get one MSE per model/cluster
        mse_model = zeros(n_clusters, 1);
        for ii = 1 : n_clusters
            % MSE check
            temp = abs(obj.local_pca{ii}.mse(:,jj)) ./ (abs(obj.local_pca{ii}.kriged_cpca_scores(:,jj)) + eps);
            mse_model(ii) = max( mean( temp, 1) ); % MSE for the ii-th model
        end
        sol = zeros(dim_out, n_clusters); % Columns are each model's prediction for jj-th prediction point 
        for ii = 1 : n_clusters
            sol(:,ii) = obj.local_pca{ii}.kcpca_predictions(:,jj); % Load prediction for model ii
        end
        % Weights inversely proportional to MSE
        w_mse = 1 ./ mse_model; 
        w_mse = w_mse ./ sum(w_mse); % Weights for the jj-th prediction
        % Weights based on distnace
        w_dist = 1 ./ (dist(:,jj) + a_tol);
        w_dist = w_dist ./ sum(w_dist);
        % Choose how to combine/get predictions
        if obj.cluster_prediction == 1
            % Weighed average of the predictions
            obj.klcpca_predictions(:,jj) = sol * (w_mse);
        elseif obj.cluster_prediction == 2
            % Use only the distance as meta-feature
            obj.klcpca_predictions(:,jj) = sol * (w_dist);
        elseif obj.cluster_prediction == 3
            % Use also the distances as meta-features
            w = (w_mse + w_dist) / 2;
            obj.klcpca_predictions(:,jj) = sol * w;
        elseif obj.cluster_prediction == 0
            % Get the closest one, instead of combining
            [~, I] = max(w_dist); % Cluste that will predict this one
            obj.klcpca_predictions(:,jj) = sol(:,I);
        end
    end
end
% They predicted the centered-scaled data
obj.klcpca_predictions = obj.uncenterUnscale(... 
                        obj.klcpca_predictions, ...
                        obj.mean_column, ...
                        obj.scaling_factors);    
end



