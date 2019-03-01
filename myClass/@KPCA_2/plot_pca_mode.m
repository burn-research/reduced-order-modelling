function [varargout] = plot_pca_mode(obj, mode_idx, var, varargin)
%% Description
% Plot the PCA mode(s)
%

%% Input
if strcmp(var, 'all')
    modes = obj.pca_modes(:,mode_idx);
    xlab = 'Variable index $[-]$';
else
    modes = obj.get_variable(mode_idx, var, 'pc', obj.is_mesh_variable);
    xlab = 'x $[cm]$';
end
[N, n_modes] = size(modes);

% Find variable's name's position
var_i = 0;
for ii = 1 : length(obj.variable_names)
    if strcmp(var, obj.variable_names{ii})
        var_i = ii;
    end
end
if var_i == 0
    x = 1:N;
else
    x = obj.mesh{var_i};
end

if iscell(obj.mesh)
    mesh_dim = size(obj.mesh{1},2);
else
    mesh_dim = size(obj.mesh,2);
end

%% Main
lgnd = {};
Markers = {'*' , '+' , 'o' , 'x', 's', 'd'};
mk = 1;
if ~obj.is_mesh_variable
    figure();
    for i = 1 : n_modes
        plot(1:N, modes(:,i), 'LineStyle','--', 'Marker',Markers{mk}); hold on;
        lgnd{end+1} = ['mode ', num2str(i)];
        % Marker
        mk = mk + 1;
        if mk > length(Markers)
            mk = 1;
        end
    end
    grid on;
    xlabel('Variable index $[-]$','Interpreter','LaTex');
    ylabel('$[-]$','Interpreter','LaTex');
    title('\textbf{PCA modes}','Interpreter','LaTex');
    return
end

if mesh_dim == 1
    figure();
    for i = 1 : n_modes
        plot(x, modes(:,i), 'LineStyle','--', 'Marker',Markers{mk}); hold on;
        lgnd{end+1} = ['mode ', num2str(i)];
        % Marker
        mk = mk + 1;
        if mk > length(Markers)
            mk = 1;
        end
    end
    grid on;
    xlabel(xlab,'Interpreter','LaTex');
    ylabel('$[-]$','Interpreter','LaTex');
    title(['\textbf{PCA modes: ', var,'}'],'Interpreter','LaTex');
    legend(lgnd);
    return
end

end

