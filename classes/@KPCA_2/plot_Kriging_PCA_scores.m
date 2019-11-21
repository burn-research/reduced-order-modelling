function plot_Kriging_PCA_scores(obj, varargin)
%% Description:
% Plot the Kriging of the PCA scores.

%% Input
% idx: PCA score to plot

idx = 1;
if ~isempty(varargin)
    idx = varargin{1};
end

is_error = false;
if length(varargin) > 1
    is_error = varargin{2};
end

% Check the dimension of the input space
if ~obj.is_mesh_variable
    dim = size(obj.training_points, 2);
else
    dim = size(obj.training_points, 2) - size(obj.mesh,2);
end

% Stop if dim > 2
if dim > 2
    error('dim > 2. Cannot plot.');
end

%% Plot
% Add label
label_string = [num2str(idx), '-th PCA score'];

% Plot routines
figure();
if ~is_error
    % Plot Kriging surfaces
    if ~isempty(obj.pca_scores_kriging_targets)
        plot_Kriging_surface(obj.training_points, obj.prediction_points, ...
            obj.pca_scores(idx,:), obj.kriged_pca_scores(idx,:), label_string, ...
            obj.pca_scores_kriging_targets(idx,:));
    else
        plot_Kriging_surface(obj.training_points, obj.prediction_points, ...
            obj.pca_scores(idx,:), obj.kriged_pca_scores(idx,:), label_string);
    end
    if ~isempty(obj.cpca_scores) && ~isempty(obj.kriged_cpca_scores) 
        figure(); label_string = [num2str(idx), '-th CPCA score'];
        plot_Kriging_surface(obj.training_points, obj.prediction_points, obj.cpca_scores(idx,:), obj.kriged_cpca_scores(idx,:), label_string);
    end
else
    % Plot regression absolute errors
    temp = obj.kriged_pca_scores(idx,:)- obj.pca_scores_kriging_targets(idx,:);
    plot_Kriging_surface(obj.training_points, obj.prediction_points, ...
        obj.pca_scores(idx,:), temp, label_string);
    title('Absolute Error [-]');
    % Plot regression relative errors
    temp = obj.kriged_pca_scores(idx,:)- obj.pca_scores_kriging_targets(idx,:);
    temp = 100 * temp ./ (abs(obj.pca_scores_kriging_targets(idx,:)) + eps);
    figure();
    plot_Kriging_surface(obj.training_points, obj.prediction_points, ...
        obj.pca_scores(idx,:), temp, label_string);
    title('Relative Error [%]');
end

end

