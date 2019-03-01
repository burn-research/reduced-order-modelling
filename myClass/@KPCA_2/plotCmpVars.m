function [varargout] = plotCmpVars(obj, varargin)
%% Descritpion
%{

%}

%% Input
A_orig = obj.original_data;

nin = length(varargin); % Number of additional inputs


% Ask to show the variables only
if nin > 0 && isa(varargin{1},'char') && strcmp(varargin{1},'show')
    pickItUp(obj.variable_names);
    return;
end


% 1st Input: Parameter
p = [];
if nin > 0
    p = varargin{1};
end


% 2nd Input: Variable
vars = obj.variable_names;
if nin > 1
    i = varargin{2};
    if ~isa(i,'char')
        vars = obj.variable_names{i};
    else
        vars = i;
    end
end
if isempty(vars) % This way you can skip this input (pass an empty vector)
    vars = obj.variable_names;
end


% 3rd Input: Avoid plotting the KLPCA
iLocal = true;
if nin > 2
    if ~isempty(varargin{3})
        iLocal = varargin{3};
    end
    if ~islogical(iLocal) || ~ismember(iLocal, [0, 1])
        error('Third input must be LOGICAL.');
    end
end


% 4th Input: Plot error
plot_error = false;
if nin > 3
    plot_error = varargin{4};
    if ~plot_error(iLocal) || ~ismember(plot_error, [0, 1])
        plot_error = false;
    end
end


%% Plot
[~, ~, y_original] = plotThisVar([], vars, A_orig, obj.map_rows, obj.xp_sorted, p, dim);

hold on;

[~, ~, y_kpca] = plotThisVar([], vars, obj.Y_sorted, obj.map_rows, obj.xp_sorted, p, dim);


xlabel('x');
if size(xq,2) == 2
    ylabel('y');
end
legend('Original', 'KPCA');
if obj.con
    legend('Original', 'KCPCA');
end


% Evaluate the errors
err_kpca = relErr(y_kpca, y_original);
err_klpca = [];


% If Local PCA has been applied
% pause(1); % Give Matlab the time to plot what it can now
if length(obj.lpca) > 1 && iLocal
    [~, ~, y_klpca] = plotThisVar(xq, vars, obj.Y_lpca, obj.map_rows, obj.xp_sorted, p, dim);
    legend('Original', 'KPCA', 'KLPCA');
    if obj.con && obj.lpca{1}.con
        legend('Original', 'KCPCA', 'KLCPCA');
    elseif obj.con && ~obj.lpca{1}.con
        legend('Original', 'KCPCA', 'KLPCA');
    end
    err_klpca = relErr(y_klpca, y_original); % Evaluate the error
end


if ~isempty(p)
    t = [vars,' - ',p_name,': [',num2str(p),']'];
    title( t );
end

% Change lines style (if 1D)
if size(obj.xq, 2) < 2
    hline = findobj(gcf, 'type', 'line');
    for i = 1 : length(hline) - 1
        set(hline(i),'LineStyle','--');
    end
end

hold off;



if plot_error
    this_y = relErr(y_kpca, y_original, true);
    figure();
    plot( this_y ); grid on;
    if length(obj.lpca) > 1 && iLocal
        hold on;
        this_y = relErr(y_kpca, y_original, true);
        plot( this_y ); 
        hold off;
    end
    xlabel(' '); xlim([1, length(y_original)]);  
    ylabel('Error [-]');
    title('Relative Error');
end

% disp([err_kpca, err_klpca]);


%% Output
if nargout > 0
    varargout{1} = err_kpca;
end
if nargout > 1
    varargout{2} = err_klpca;
end



end



