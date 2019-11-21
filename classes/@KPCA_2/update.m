function update(obj, varargin)
 
% Load object properties
c = properties(obj);


% NO INPUT CASE: Update the property whose name is in STR
if isempty(varargin)
    obj.updateALL();
    obj.str = [];
    return;
end




% ONE OR TWO INPUTS CASE: 1st input is the value, 2st ipunt (optional) the
% property to update with that value

v = varargin{1};
if nargin > 1
    % Load user-provided string
    s = varargin{2};    
else
    % Ask the user to input the string
    i = pickItUp(c, '\n Which property to will trigger the updates? \n'); 
    s = c{i};
end

% Check the string is an actual property
i = strcmp(s,c);

% If so, update that property with the user-provided value V
if any(i)
    prompt = ['obj.',s,' = v;'];
    eval(prompt);
else
    error('No property with name %s.',s);
end

% Update all dependet variables
obj.str = s;
obj.updateALL();
obj.str = [];
        



end




    