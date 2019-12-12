function plot_sampled_space(varargin)
%% Description

%% Main
% Varargin
n_args = length(varargin);
leg = {}; % Legend 
M = {'o', '*', '+', 's'}; mk = 1;
C = {'k', 'r', 'g', 'b', 'm'}; cl = 1;
figure();
if n_args > 0
    for ii = 1 : n_args
        if ~isempty(varargin{ii})
            inner_plot(varargin{ii}, M{mk}, C{cl});
            leg{end+1} = num2str(ii);
            mk = mk + 1;
            cl = cl + 1;
            if mk > length(M)
                mk = 1;
            end
            if cl > length(C)
                cl = 1;
            end
        end
    end
    % Legend and title
    legend(leg);
    title('Sampled space');
    hold off;
    return;
elseif n_args == 0
    fprintf('No input provided.\n');
    return;
end

end

function inner_plot(setofpoints, mk, cl)
% Check the points are there 
if isempty(setofpoints) 
    fprintf('No points are available.\n');
    return;
end
% Check the dimension of the input space
dim = size(setofpoints, 2);
if dim > 3
    fprintf('Dimension of input space is > 3.\nNo plot possible.\n');
    return;    
end
% Input space with dimension = 1
if dim == 1
    y_vec = setofpoints * 0 + 1;
    plot(setofpoints, y_vec, 'Marker',mk, 'Color',cl, 'LineStyle','none');
    hold on; grid on;
end
% Input space with dimension = 2
if dim == 2
    plot(setofpoints(:,1), setofpoints(:,2), 'Marker',mk, 'Color',cl, 'LineStyle','none');
    hold on; grid on;
end
% Input space with dimension = 3
if dim == 3 
    plot3(setofpoints(:,1), setofpoints(:,2), setofpoints(:,3), 'Marker',mk, 'Color',cl, 'LineStyle','none');
    hold on; grid on;
end
end



