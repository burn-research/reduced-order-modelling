function plot_compare_variable(obj, var, point, varargin)
%% Input

n_args = length(varargin);

% Input: var
if ~isa(var,'char')
    var = obj.variable_names{var};
end

%% Main

% Initialize legend and title
lgnd = {'Data'};
ttl = ['\textbf{',var,'}'];

% Get vectors to plot
y_data = obj.get_variable(point, var, 'data', obj.is_mesh_variable);
y_direct = obj.get_variable(point, var, 'direct', obj.is_mesh_variable);
y_pca = obj.get_variable(point, var, 'pca', obj.is_mesh_variable);
y_cpca = obj.get_variable(point, var, 'cpca', obj.is_mesh_variable);
y_lpca = obj.get_variable(point, var, 'lpca', obj.is_mesh_variable);
y_lcpca = obj.get_variable(point, var, 'lcpca', obj.is_mesh_variable);

% Get x-axis data
if size(obj.mesh, 2) == 1
    x = obj.mesh;
else
    x = 1:length(y_data); 
end

% Plot
figure();
s1 = '-'; s2 = '--';
if n_args > 0
    s1 = varargin{1};
end
if n_args > 1
    s2 = varargin{2};
end
plot(x, y_data, 'LineStyle', s1, 'Color','k'); 
hold on;
if ~isempty(y_direct)
    plot(x, y_direct, s2);
    lgnd{end+1} = 'Direct'; % Append to the legend
end
if ~isempty(y_pca)
    plot(x, y_pca, s2);
    lgnd{end+1} = 'PCA'; % Append to the legend
end
if ~isempty(y_cpca)
    plot(x, y_cpca, s2);
    lgnd{end+1} = 'CPCA'; % Append to the legend
end
if ~isempty(y_lpca)
    plot(x, y_lpca, s2);
    lgnd{end+1} = 'LPCA'; % Append to the legend
end
if ~isempty(y_lcpca)
    plot(x, y_lcpca, s2);
    lgnd{end+1} = 'LCPCA'; % Append to the legend
end

% Plot settings
xlabel('Point index $[-]$','Interpreter','LaTex');
ylabel('Value $[-]$','Interpreter','LaTex');
grid on;
legend(lgnd);
title(ttl,'Interpreter','LaTex');
hold off;


end



