function plot_Influence(points, Infl, varargin)

% Check the dimension of the input space
dim = size(points, 2);
if dim > 3
    fprintf('Dimension of input space is > 3.\nNo plot possible.\n');
    return;    
end

% Number of points 
n_points = size(points, 1);

% Symbols
s1 = 'ko';
s2 = 'bx';

% Start new figure
figure();

% Input space with dimension = 1
if dim == 1
    semilogy(points, 100 * Infl, '--*k');
    grid on;
    xlabel('Parameter value $[-]$','Interpreter','LaTex');
    ylabel('Relative Influece $[\%]$','Interpreter','LaTex');
end

% Input space with dimension = 2
if dim == 2
    plot(points(:,1), points(:,2), s1);
    hold on; grid on;
    plot(points(:,1), points(:,2), s2);
    hold on; grid on;
    % Plot errors
    for ii = 1 : n_points
        x = points(ii,1);
        y = points(ii,2);
        txt = [' ', num2str(100 * Infl(ii), 3), '%'];
        text(x, y, txt);
    end
end

% Input space with dimension = 3
if dim == 3 
    plot3(points(:,1), points(:,2), points(:,3), s1);
    hold on; grid on;
    plot3(points(:,1), points(:,2), points(:,3), s2);
    hold on; grid on;
    % Plot errors
    for ii = 1 : n_points
        x = points(ii,1);
        y = points(ii,2);
        z = points(ii,3);
        txt = [' ', num2str(100 * Infl(ii), 3), '%'];
        text(x, y, z, txt);
    end
end

% Legend and title
tit = ['\textbf{Relative Influence}'];
title(tit, 'Interpreter', 'LaTex');
hold off;

end


