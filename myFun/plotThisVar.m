function [varargout] = plotThisVar(x, v, Y, map, xp, varargin)
%% Descritpion
%{
Y is a multivariate data matrix. We are interested in plotting the variable
V for the case P, against the data X.
MAP and XP are needed to get the right row(s) and column(s) from Y.
%}


%% Input
[rows, cols] = size(Y);
if size(xp,1) ~= cols
    error('Number of columns of Y must be equal to rows of XP.');
end


% Varargin
p = [];
nin = length(varargin);
if nin > 0
    p = varargin{1};
end

dim = 0;
if nin > 1
    dim = varargin{2};
end
if isempty(dim)
    dim = 0;
end


% v
if ~isa(v,'char') && isa(v,'cell')
    prompt = '\nWhich variable?\n';
    vars = v;

    i = pickItUp(vars,prompt);
    v = vars{i};
end



%% Main
[this_rec, i] = getThisVar(Y, map, v);

% Get the case
j = [];
if isempty(p)
    prompt = '\nFor which case?\n'; % Use personalizes prompt
    p = xp(:,dim+1:end); % Load the list of parameter points
    j = pickItUp(p, prompt); % Get the row/point you need
    p = p(j,:); % Load that row/point
end


% If also specific columns have to be selected
temp = xp(:,dim+1:end);
I = ismember(temp, p, 'rows');
y = this_rec(:, I);



%% Varargout
if nargout > 0
    varargout{1} = i;
    if nargout > 1
        varargout{2} = j;
    end
    if nargout > 2
        varargout{3} = y;
    end
end


%% Plot
setFigOpts;

if ~isempty(x)
    if size(x,2) == 2
      % 2D Mesh
        plot3(x(:,1), x(:,2), y, '.');
%         surfaceFrom3Vectors(x(:,1), x(:,2), y);
      %%%
    elseif size(x,2) == 1
      % 1D Mesh
        plot(x,y);
      %%%
    else
        warning('2nd dimension of the mesh was neither 1 nor 2.');
        plot(x,y);
    end
else
    plot(y);
end
grid on;
title(v);


end

