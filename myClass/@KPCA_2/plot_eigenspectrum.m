function plot_eigenspectrum(obj, varargin)

% Default values
I = 1:length(obj.pca_eigenvalues);

% Store the value of scaling_criterion upon entering this routine
sc_entry = obj.scaling_criterion;

% Inputs
n_args = length(varargin);
if n_args > 0 && ~isempty(varargin{1})
    % Eigenvalues to plot
    I = varargin{1};
end
all_sc = false;
if n_args > 1
    if ischar(varargin{2}) && (strcmp(varargin{2}, 'all_sc') || strcmp(varargin{2}, 'sc_all'))
        all_sc = true;
    else 
        all_sc = varargin{2};
    end
end

if ~all_sc
    % Eigenvalues to plot
    plt = obj.pca_eigenvalues(I);
    % Plot
    figure();
    semilogy(plt, 'LineStyle','--', 'Color','k', 'Marker','*'); 
else
    figure();
    Markers = {'*', 'x', '+', 'd', 'p', 'o'};
    mk = 1;
    leg = {};
    for sc = 1 : 6
        obj.scaling_criterion = sc;
        obj.runPca();
        % Eigenvalues to plot
        plt = obj.pca_eigenvalues(I); hold on;
        % Plot
        semilogy(plt, 'LineStyle','--', 'Marker',Markers{mk});
        leg{end+1} = ['sc: ', num2str(sc)];
        % Marker update
        mk = mk + 1;
        if mk > length(Markers) 
            mk = 1;
        end
    end
    legend(leg); % Legend
    % Reset initial status
    obj.scaling_criterion = sc_entry;
    obj.runPca();
end
set(gca, 'YScale', 'log');
grid on;
% Axis labels
xlabel('Index of the eigenvalues','Interpreter','LaTex');
ylabel('Value [-]','Interpreter','LaTex');
% Title 
title('\textbf{Eigenvalue spectrum}','Interpreter','LaTex');
hold off;


end











