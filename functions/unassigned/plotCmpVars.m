function plotCmpVars(x, v, Y, Y_orig, row_map, col_map, varargin)
%% Description
%{ 
Given a dataset Y_orig that we tried to predict with the KPCA approach, now
we want to compare graphically the predicted variables' profiles with the
actual ones.
This is basically like running the plotThisVar function twice, with an HOLD
ON.

INPUT
x: mesh
v: variable or list of variables (from which to pick one up)
Y: KPCA dataset
Y_orig: actual dataset
row_map: map for the rows
col_map: map for the cols

OUTPUT
None.

%}

%% Input

% Varargin
p = [];
nin = length(varargin);
if nin > 0
    p = varargin{1};
end

% Get the case
if isempty(p)
    prompt = '\nFor which case?\n';
    p = unique(col_map(:, size(col_map,2) ));

    j = pickItUp(p,prompt);
    p = p(j);
end



%% Plot
i = plotThisVar(x, v, Y_orig, row_map, col_map, p);
hold on;
plotThisVar(x, v{i}, Y, row_map, col_map, p);

legend('Original','KPCA');




end


