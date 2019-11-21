function [varargout] = plot_clustered_variable_space(obj, char_var, varargin)
%% Description:
% Plot how the mesh-points have been clustered for a specific variable.
% Input
    % (1) Physical variable to be analysed, or input space point.
%

%% Input
% Dummy check
if obj.is_mesh_variable && obj.clustering_dimension
    error('[plot_clustered_variable_space] Not possible: kpca.is_mesh_variable && kpca.clustering_dimension');
end

n_args = length(varargin);

% Input: var
if ~isa(char_var, 'char')
    char_var = obj.variable_names{char_var};
end

% Input: boolean_var
x = obj.training_points;

% Input: par (e.g. eq. ratio, speed, etc.)
if n_args > 0
    par = n_args;
else
    par = 1;
end

%% Main
m_size = size(obj.mesh, 2);

% Get vectors to plot
y = vectorize(obj.get_variable(x, char_var, 'data', obj.is_mesh_variable));
x = x(:, m_size+par); % Get rid of the mesh values and get the desired par

% Plot
M = {'o', '+', '*', '.', 'x'};
C = {'k', 'r', 'b', 'g', 'y', 'm', 'c'};
n_clusters = max(obj.local_idx);
figure(); hold on;
leg = {};
k = 1; % Markers iterator
c = 1; % Colors iterator
for ii = 1 : n_clusters
    I = ( obj.local_idx == ii );
    plot(x(I,:), y(I), 'Marker', M{k}, 'Color', C{c}, 'LineStyle','none');
    leg{end+1} = ['Cluster', num2str(ii)];
    k = k + 1;
    c = c + 1;
    if k == length(M)
        k = 1; % Reset markers iterator
    end
    if c == length(C)
        c = 1; % Reset colors iterator
    end
end
legend(leg);
% xlabel('$ \mathbf{Input \ point \ [-]} $', 'Interpreter','LaTex');
% ylabel(['$ \mathbf{', char_var,'} $'], 'Interpreter','LaTex')
title([char_var]);
grid on;
hold off;

end

% Local functions
function y = vectorize(X)

y = zeros(size(X,1) * size(X,2), 1);
for ii = 1 : size(X, 2)
    y(1+(ii-1)*size(X,1):(ii*size(X,1))) = X(:,ii);
end

end


