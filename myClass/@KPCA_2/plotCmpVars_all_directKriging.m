function [varargout] = plotCmpVars_all_directKriging(obj, varargin)
%% Description
%{

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


%% Main

if k == 0
    % Training points
    err_directKriging = cell(size(obj.xp,1),1);     
    for i = 1 : size(obj.xp,1)
        figure();
        err_directKriging{i} = obj.plotCmpVars_directKriging(obj.xp(i,:), input_string);
        hold on; obj.plotCmpVars(obj.xp(i,:), input_string);
        children = get(gca, 'children'); delete(children(3));
        legend('Original','Direct Kriging','KPCA','KLPCA');
        pause(this_time);
    end
elseif k == 1
    % Prediction points
    err_directKriging = cell(size(obj.xp_kriged,1),1);     
    for i = 1 : size(obj.xp_kriged,1)
        figure();
        err_directKriging{i} = obj.plotCmpVars_directKriging(obj.xp_kriged(i,:), input_string);
        hold on; obj.plotCmpVars(obj.xp_kriged(i,:), input_string);
        children = get(gca, 'children'); delete(children(3));
        legend('Original','Direct Kriging','KPCA','KLPCA');
        pause(this_time);
    end
end


%% Output

if nargout > 0 
    % Direct Kriging errors
    varargout{1} = zeros(length(err_directKriging),1);
    for i = 1 : length(err_directKriging)
        varargout{1}(i,1) = mean(err_directKriging{i});
    end
end


end









