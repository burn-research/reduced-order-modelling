function plotLpcaIdx(obj, varargin)
%% Descritpion
%{

%}

%% Input
nin = length(varargin);


% Variable
vars = obj.vars;
if nin > 0
    i = varargin{1};
    if ~isa(i,'char')
        vars = obj.vars{i};
    else
        vars = i;
    end
else
    i = pickItUp(obj.vars, 'Choose the variable to plot:\n');
    vars = obj.vars{i};
end
% This way you can skip this input (pass an empty vector)
if isempty(vars) 
    vars = obj.vars;
end


% Get the correct grid
xq = obj.xq;
if isa(xq,'cell')
    xq = [];
end



%% Main
[selected_data, ~, I] = getThisVar(obj.Y, obj.map_rows, vars);


%% Plot
y2plot = obj.idx(I); 
plot(y2plot, '.'); grid on;
xlabel('Grid points');  xlim([1 length(y2plot)]);
ylabel('Cluster index');
title(vars);

end


