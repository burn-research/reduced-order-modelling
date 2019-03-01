function plot_parity(obj, char_var, boolean_var, varargin)
%% Description:
% Plot the observed vs recovered data.
% 
% Inputs
% (1) char_var: Physical variable to be analysed.
% (2) boolean_var: (TRUE/FALSE) plot the reconstruction or prediction case.
% (2) numerical_var: User wants to plot for this point.
%

%% Input
% Input: var
if ~isa(char_var, 'char')
    char_var = obj.variable_names{char_var};
end

% Input: boolean_var
if isempty(boolean_var) 
    points = obj.prediction_points;
elseif islogical(boolean_var) && boolean_var
    points = obj.prediction_points;
elseif islogical(boolean_var) && ~boolean_var
    points = obj.training_points;
elseif isnumeric(boolean_var)
    points = boolean_var;
    boolean_var = ismember(points, obj.prediction_points, 'rows');
else 
    error('Second input must be either LOGICAL or NUMERIC.');
end

% Input: smooth
smoothing = false;
if ~isempty(varargin) && islogical(varargin{1})
    smoothing = varargin{1};
end

% Input: n_smooth
n_smooth = 1;
if length(varargin) > 1
    n_smooth = varargin{2};
    if n_smooth < 1
        smoothing = false;
    end
end


%% Main

% Get vectors to plot
y_data = vectorize(obj.get_variable(points, char_var, 'data', obj.is_mesh_variable));
y_pca = vectorize(obj.get_variable(points, char_var, 'pca', obj.is_mesh_variable));
y_cpca = vectorize(obj.get_variable(points, char_var, 'cpca', obj.is_mesh_variable));
y_lpca = vectorize(obj.get_variable(points, char_var, 'lpca', obj.is_mesh_variable));
y_lcpca = vectorize(obj.get_variable(points, char_var, 'lcpca', obj.is_mesh_variable));
y_direct = vectorize(obj.get_variable(points, char_var, 'direct', obj.is_mesh_variable));

% Smoothing
if smoothing
    for i = 1 : n_smooth
        y_pca = smooth(y_data, y_pca);
        y_cpca = smooth(y_data, y_cpca);
        y_lpca = smooth(y_data, y_lpca);
        y_lcpca = smooth(y_data, y_lcpca);
        y_direct = smooth(y_data, y_direct);
    end
end

% Normalizing factor
norm_fac = 1; % abs(max(y_data));

% Plot
Markers = {'o' , '*' , 'x' , '+', 's'};
s = ' ';
leg = {};
hf = figure(); hold on; grid on; 
if ~isempty(y_pca)
    plot(y_data/norm_fac, y_pca/norm_fac, 'MarkerSize',8, 'Marker',Markers{1}, 'Color',[.65 .65 .65], 'LineStyle','none');
    err = get_err(y_data, y_pca);
    leg = {['pca ', s]};
    err_bar = [min([y_data(:); y_pca(:)]), max([y_data(:); y_pca(:)])];
end
if ~isempty(y_cpca)
    hold on;
    plot(y_data/norm_fac, y_cpca/norm_fac, 'MarkerSize',8, 'Marker',Markers{2}, 'Color','b', 'LineStyle','none'); 
    err = get_err(y_data, y_cpca);
    leg{end+1} = ['cpca ', s];
    err_bar = [err_bar; min(y_cpca(:)), max(y_cpca(:))];
end
if ~isempty(y_lpca)
    hold on;
    plot(y_data/norm_fac, y_lpca/norm_fac, 'MarkerSize',8, 'Marker',Markers{3}, 'Color','k', 'LineStyle','none');
    err = get_err(y_data, y_lpca);
    leg{end+1} = ['lpca ', s];
end
if ~isempty(y_lcpca)
    hold on;
    plot(y_data/norm_fac, y_lcpca/norm_fac, 'MarkerSize',8, 'Marker',Markers{4}, 'Color','m', 'LineStyle','none'); 
    err = get_err(y_data, y_lcpca);
    leg{end+1} = ['lcpca ', s];
    err_bar = [err_bar; min(y_lcpca(:)), max(y_lcpca(:))];
end
if ~isempty(y_direct)
    hold on;
    plot(y_data/norm_fac, y_direct/norm_fac, 'MarkerSize',8, 'Marker',Markers{5}, 'Color','g', 'LineStyle','none'); 
    err = get_err(y_data, y_direct);
    leg{end+1} = ['direct ', s];
    err_bar = [err_bar; min(y_direct(:)), max(y_direct(:))];
end
% Plot x = y line
% plot([min(y_data) max(y_data)]/norm_fac, [min(y_data) max(y_data)]/norm_fac, 'y--'); 
% Plot error lines
val1 = min(y_data);
val2 = max(y_data);
a = .95;
plot([val1 val2], [val1*a val2*a], '--r', 'LineWidth', 2);
a = 1.05;
plot([val1 val2], [val1*a val2*a], '--r', 'LineWidth', 2);
a = .9;
plot([val1 val2], [val1*a val2*a], '-.', 'LineWidth', 2, 'Color',[.93 .69 .13]);
a = 1.1;
plot([val1 val2], [val1*a val2*a], '-.', 'LineWidth', 2, 'Color',[.93 .69 .13]);
% Title and labels
% tit = ['\textbf{',char_var, '}'];
% title(tit, 'Interpreter','LaTex');
if norm_fac ~= 1
    xlabel('Normalized observed value ', 'Interpreter','LaTex');
    if ~boolean_var
        ylabel('Normalized recovered value ', 'Interpreter','LaTex');
    else
        ylabel('Normalized predicted value ', 'Interpreter','LaTex');
    end
else
    xlabel(['Observed value ', char_var], 'Interpreter','LaTex');
    if ~boolean_var
        ylabel(['Recovered value ', char_var], 'Interpreter','LaTex');
    else
        ylabel(['Predicted value ', char_var], 'Interpreter','LaTex');
    end
end
% Legend
legend(leg);
%figure_name = [char_var, '_', date];
%figure_extension = 'eps';
%saveas(gcf, figure_name, figure_extension);
%print(hf, '-dpdf', [figure_name,'.pdf'], '-opengl');
set(findall(gcf,'-property','FontSize'),'FontSize',25);
hold off;

end

% Local functions
function y = vectorize(X)

y = zeros(size(X,1) * size(X,2), 1);
for ii = 1 : size(X, 2)
    y(1+(ii-1)*size(X,1):(ii*size(X,1))) = X(:,ii);
end

end

function err = get_err(y1, y2)

z = sqrt( (y2 - y1).^2 );
err = mean( z ) / mean(y1);

end





