function [varargout] = massFracCheck(obj, varargin)


%% Input
nin = length(varargin);





%% Main

% Get the data and the variables related to chemical species
this = obj.Y_sorted(2*obj.nq + 1:end, :);
this_lpca = obj.Y_lpca(2*obj.nq + 1:end, :);
theseVars = getTheVars(obj.vars);


% Change the shape of the data matrix
this = getStateVars( this, obj.xq, theseVars, [] );
this_lpca = getStateVars( this_lpca, obj.xq, theseVars, [] );


% Evaluate the total mass fraction for each case in every grid point
t = size(this,2) / length(obj.xq);
cases = cell(t,1); cases_lpca = cases; 
totalMassFraction = cases; totalMassFraction_LPCA = cases;
for i = 1 : t
    % For KPCA
    cases{i} = this(:, 1+(i-1)*obj.nq:i*obj.nq);
    totalMassFraction{i} = sum(cases{i}, 1);

    % For KLPCA
    cases_lpca{i} = this_lpca(:, 1+(i-1)*obj.nq:i*obj.nq);
    totalMassFraction_LPCA{i} = sum( cases_lpca{i}, 1);
end



%% Output
if nargout > 0
    varargout{1} = totalMassFraction;
end
if nargout > 1
    varargout{2} = totalMassFraction_LPCA;
end



end



function v = getTheVars(this)

n = 2;

for i = 1 : length(this) - n
    v{i,1} = this{i+n};
end


end









