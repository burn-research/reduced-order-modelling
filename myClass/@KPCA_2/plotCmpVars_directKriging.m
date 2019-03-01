function [varargout] = plotCmpVars_directKriging(obj, varargin)
%% Descritpion
%{

%}

%% Strings for the plot
stringPlot = cell(3, 1);
stringPlot{1} = 'Original';
stringPlot{2} = 'Direct Kriging';
stringPlot{3} = 'KPCA';
stringPlot{4} = 'KCPCA';


%% Input
A_orig = obj.Y_orig;

nin = length(varargin); % Number of additional inputs

% Ask to show the variables only
if nin > 0 && isa(varargin{1},'char') && strcmp(varargin{1},'show')
    pickItUp(obj.vars);
    return;
end

% 1st Input: Parameter
p = [];
if nin > 0
    p = varargin{1};
end

% 2nd Input: Variable
vars = obj.vars;
if nin > 1
    i = varargin{2};
    if ~isa(i,'char')
        vars = obj.vars{i};
    else
        vars = i;
    end
end
if isempty(vars) % This way you can skip this input (pass an empty vector)
    vars = obj.vars;
end

% Get the correct grid
xq = obj.xq;
if isa(xq,'cell')
    xq = [];
end

% Mesh size
dim = size(obj.xq, 2);
if ~strcmp(obj.x_status,'parameter')
    dim = 0;
end

% Parameter name
p_name = 'Point';
if ~isempty(obj.xp_names) 
    p_name = obj.xp_names{1};
end


%% Plot
[~, ~, y_original] = plotThisVar(xq, vars, A_orig, obj.map_rows, obj.xp_sorted, p, dim);

hold on;

[~, ~, y_directKriging] = plotThisVar(xq, vars, obj.Y_k_d_sorted, obj.map_rows, obj.xp_sorted, p, dim);


xlabel('x');
if size(xq,2) == 2
    ylabel('y');
end
legend(stringPlot{1}, stringPlot{2});


% Evaluate the errors
err_directKriging = relErr(y_directKriging, y_original);


if ~isempty(p)
    t = [vars,'    - ',p_name,': [',num2str(p),']'];
    title( t );
end


% Change lines style
hline = findobj(gcf, 'type', 'line');
for i = 1 : length(hline) - 1
    set(hline(i),'LineStyle','--');
end


hold off;


%% Output
if nargout > 0
    varargout{1} = err_directKriging;
end
if nargout > 1
    varargout{2} = err_klpca;
end


end


