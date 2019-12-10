function i = pickItUp(list, varargin)
%% Description
% Ask the user which variable/case he wants to pick up.

%% Build question

% Initialize the question
prompt0 = '\n Pick it up: \n'; prompt=[];

% Varargin
if ~isempty(varargin)
    prompt0 = varargin{1};
end

% Size of LIST
t = size(list,1);

for i = 1 : t    
    % Class of LIST
    if isa(list,'cell')
        this = list{i};
    elseif isa(list,'double')
        this = num2str(list(i));
    end
    
    % Append
    prompt = [prompt,'[',num2str(i),']','-',this,'   '];
    if (i/8) == round(i/8)
        prompt = [prompt,'\n'];
    end
    
end

%% Ask question
i = input([prompt0, prompt,'\n \n']);



end

