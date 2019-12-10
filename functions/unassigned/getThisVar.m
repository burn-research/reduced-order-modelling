function [selected_data, i, I] = getThisVar(A, map, v, varargin)
%% Description
%{
The variables is searched for along the first dimension of A.

INPUT
map: cell array
var: string

%}


%% Input

% Size
[rows, cols] = size(A);
[map_rows, map_cols] = size(map);


% Size check
if rows ~= map_rows
    error('Sizes dont match.');
end


% Get the column of MAP that contains STRINGs
c = 1;
while map_cols > 1 && ~isa(map{1,c},'char')
    c = c + 1;
end


i = [];
if ~isa(v,'char') && isa(v,'cell')
    vars = v;
    
    % Build question
    prompt0 = ['Choose the variable to get: \n']; prompt = [];
    
    t = size(vars,1); 
    for i = 1 : t
        prompt = [prompt,'[',num2str(i),']','-',vars{i},'   '];
        if (i/8) == round(i/8);   prompt = [prompt,'\n'];   end
    end
    
    % Ask question
    i = input([prompt0, prompt,'\n \n']);
    v = vars{i};
end



%% Main
I = false(map_rows, 1);
for l = 1 : map_rows
    I(l) = strcmp(v, map{l,c});
end



%% Output
selected_data = A(I,:);




end

