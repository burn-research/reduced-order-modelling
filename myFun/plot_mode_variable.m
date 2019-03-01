function plot_mode_variable(mode, mesh, variables)
%% Description
% As POD for spatial fields involves modes that account for more than one
% spatial fields (thus, more than one variable), this function is used to
% differentiate the weights on one mode corresponding to the different
% spatial fields.

%% Input
% Check
if (size(mode,1) ~= 1) && (size(mode,2) ~= 1)
    error('First input argument MODE must be a vector.');
end
% Get sizes and check inputs
mode = mode(:); % Make it column-vector
ntot = length(mode);
nx = size(mesh,1);
nvar = ntot / nx;
if nvar ~= round(nvar)
    error('Number of fields must be an integer.');
end
% Optional input (variables' names)
if nargin < 3
    variables = cell(nvar, 1);
    for ii = 1 : nvar
        variables{ii} = num2str(ii);
    end
else
    if (nvar ~= length(variables)) || ~iscell(variables)
        error('Input argument VARIABLES has wrong size or type.');
    end
end

%% Main
rgb = rand(nvar, 3);
p = cell(nvar,1);
figure(); hold on;
for ii = 1 : nvar
    x1 = 1 + (ii - 1) * nx;
    x2 = ii * nx;
    x = x1:x2;
    y = mode(x1:x2);
    p{ii} = plot(x, y, '.');
%     set(p{ii}, 'FaceColor', rgb(ii,:));
end
xlim([1 ntot])
ylim([min(mode) max(mode)]);
legend(variables);
hold off;

end







