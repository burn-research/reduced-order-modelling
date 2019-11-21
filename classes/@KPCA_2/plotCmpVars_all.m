function [varargout] = plotCmpVars_all(obj, varargin)
%% Description
%{
Ask the plotCmpVars function to plot for all cases
%}


%% Input
nin = length(varargin);


% First input
k = [];
if nin > 0
    k = varargin{1};
end
if isempty(k)
    k = input('\nTraining or Prediction points? [0/1] \n');
    if k ~= 0 && k ~= 1
        warning('Wrong input. Set to 0.'); 
        k = 0;
    end
end


% Second input
input_string = obj.vars{1};
if nin > 1
    input_string = varargin{2};
end
if ~isa(input_string,'char')
    warining('Second input must be a string. Perhaps see/provide property obj.vars. Set to obj.vars{1}.');
    input_string = obj.vars{1};
end


% Third input
this_time = 1;
if nin > 2
    this_time = varargin{3};
end


% Fourth input
worstCase = false;
if nin > 3
    worstCase = varargin{4};
    if ~islogical(worstCase) && ~ismember(worstCase,[0 1])
        worstCase = false;
    end
    if isa(worstCase,'char') && strcmp(worstCase,'worst case')
        worstCase = true;
    end
end



%% Main

if k == 0
    % Training points
    err_kpca = cell(size(obj.xp,1),1);     err_klpca = err_kpca;
    for i = 1 : size(obj.xp,1)
        figure();
        [err_kpca{i}, err_klpca{i}] = obj.plotCmpVars(obj.xp(i,:), input_string);
        pause(this_time);
    end
elseif k == 1
    % Prediction points
    err_kpca = cell(size(obj.xp_kriged,1),1);     err_klpca = err_kpca;
    for i = 1 : size(obj.xp_kriged,1)
        figure();
        [err_kpca{i}, err_klpca{i}] = obj.plotCmpVars(obj.xp_kriged(i,:), input_string);
        pause(this_time);
    end
end


%% Output

if nargout > 0 
    % KPCA errors
    varargout{1} = zeros(length(err_kpca),1);
    for i = 1 : length(err_kpca)
        varargout{1}(i,1) = mean(err_kpca{i});
    end
end

if nargout > 1
    % KLPCA errors
    varargout{2} = zeros(length(err_klpca),1);
    for i = 1 : length(err_klpca)
        varargout{2}(i,1) = mean(err_klpca{i});
    end
end


end


function [varargout] = saveit(h, iSave, obj, input_string, p, varargin)

if ~iSave
    return;
end

cn = num2str(obj.num_clusts);
k = num2str(obj.k);
p = num2str(p);

filename = [input_string,'_y_cn',cn,'_k',k,'_p',p];


savefig(h, filename);


end






