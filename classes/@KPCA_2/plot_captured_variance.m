function plot_captured_variance(obj, varargin)

n_args = length(varargin);
do_all = false;
if n_args > 0
    do_all = varargin{1};
end

if do_all 
    sc_entry = obj.scaling_criterion;
    figure(); hold on;
    Markers = {'*', 'x', '+', 'd', 'p', 'o'};
    mk = 1;
    leg = {};
    for ii = 1 : 6
        obj.scaling_criterion = ii;
        obj.runPca();
        % Get the cumulative energy
        y = [0; cumsum(obj.pca_eigenvalues)/sum(obj.pca_eigenvalues)];
        x = 0:length(y)-1;
        % Plot
        plot(x, 100 * y, 'LineStyle','--', 'Marker',Markers{mk});
        hold on;
        leg{end+1} = ['sc: ', num2str(ii)];
        % Marker update
        mk = mk + 1;
        if mk > length(Markers) 
            mk = 1;
        end
    end
    legend(leg); % Legend
    obj.scaling_criterion = sc_entry;
    obj.runPca();
else
    figure();
    % Get the cumulative energy
    y = [0; cumsum(obj.pca_eigenvalues)/sum(obj.pca_eigenvalues)];
    x = 0:length(y)-1;
    % Plot
    plot(x, 100 * y, 'LineStyle','--', 'Color','k', 'Marker','*');
end
% PCA energy
obj.pca_energy = y(obj.pca_approximation_order+1);
grid on;
% Axis scale
set(gca, 'yscale', 'log');
% Axis labels
xlabel('Number of PCs','Interpreter','LaTex');
ylabel('[$\%$]','Interpreter','LaTex');
% Axis limits
xlim([0, length(y)]);
ylim([0, 100.5]);
% Title 
title('\textbf{Captured Variance}','Interpreter','LaTex');
hold off;

end
