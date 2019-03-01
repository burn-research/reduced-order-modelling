function plot_parity_utility(varargin)
%% Description
% INPUT
%   {1} Observed data
%   {>1} Predicted/Recovered data

%% Input
n = length(varargin);
y = cell(n, 1);

%% Main
% Get vectors to plot
for ii = 1 : n
    y{ii} = vectorize(varargin{ii});
end
% Normalizing factor
norm_fac = 1; % abs(max(y{1}));
% Plot
Markers = {'o' , '*' , 'x' , '+', 's', 'd'};
Colors = {'k', [.65 .65 .65], 'b', 'm', 'g'};
leg = {};
hf = figure(); hold on; grid on; 
kk = 1;
cc = 1;
for jj = 2 : n
    plot(y{1}/norm_fac, y{jj}/norm_fac, 'MarkerSize',8, 'Marker',Markers{kk}, 'Color',Colors{cc}, 'LineStyle','none');
    leg = {['Model ', num2str(jj-1)]};
    % Markers
    kk = kk + 1;
    if kk > length(Markers)
        kk = 1;
    end
    % Colors
    cc = cc + 1;
    if cc > length(Colors)
        cc = 1;
    end
end
% Plot x = y line
% plot([min(y{1}) max(y{1})]/norm_fac, [min(y{1}) max(y{1})]/norm_fac, 'y--'); 
% Plot error lines
val1 = min(y{1});
val2 = max(y{1});
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
    xlabel('Observed value ', 'Interpreter','LaTex');
    ylabel('Predicted value ', 'Interpreter','LaTex');
end
% Legend
% legend(leg);
%figure_name = [char_var, '_', date];
%figure_extension = 'eps';
%saveas(gcf, figure_name, figure_extension);
%print(hf, '-dpdf', [figure_name,'.pdf'], '-opengl');
set(findall(hf,'-property','FontSize'),'FontSize',25);
hold off;
end
% Local functions
function y = vectorize(X)
y = zeros(size(X,1) * size(X,2), 1);
for ii = 1 : size(X, 2)
    y(1+(ii-1)*size(X,1):(ii*size(X,1))) = X(:,ii);
end
end