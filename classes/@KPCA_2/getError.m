function [varargout] = getError(obj, varargin)
%% Description
% First input: reference data
% Second input: approximated data
% Third input: corresponding input-space points (training or prediction
%              points, for example).
% The fourth and fifth inputs are optional.
%

%% Input
n_args = length(varargin); % Number of inputs
% Stop if there are not enough inputs
if n_args < 2
    fprintf('This function needs 2 inputs at least, mesh and variables.\n');
    varargout{1} = [];
    varargout{2} = [];
    return
end
% Mesh
if n_args < 3 || isempty(varargin{3})
    mesh = obj.mesh;
else
    mesh = varargin{3};
end
% Variable names (cell-array of string)
if n_args < 4 || isempty(varargin{4})
    vars = obj.variable_names;
else
    vars = varargin{4};
end

%% Evaluate errors
% The routine will access the is_mesh_variable property of the class in
% order to understand how to operate.
a_tol = 1e-8; % Tolerance
% Scaling criterion
if iscell(vars) && length(vars) > 1
    sf_C = 2;
else
    sf_C = 6;
end
if obj.is_mesh_variable
    % Difference matrix
    Delta = getError(varargin{1}, varargin{2}, 2, 0);
    % Scaling factors
    [~, ~, sf] = center_scale_data(varargin{1}, sf_C);
    % Scale the difference matrix
    Delta = Delta .* repmat(sf, 1, size(varargin{1},2)); % Multiply
    Delta = Delta ./ (repmat(sf.^2, 1, size(varargin{1},2)) + a_tol); % Divide
    % Rearrange difference matrix
    n_samples = 1:size(varargin{1},2); n_samples = n_samples(:);
    Q_delta = getStateVars(Delta, mesh, vars, n_samples);
    % First kind errors
    varargout{1} = mean(Delta, 1);
    % Second kind errors
    varargout{2} = mean(Q_delta, 2);
%     % Rearrange matrices
%     n_samples = 1:size(varargin{1},2); n_samples = n_samples(:);
%     Q_ref = getStateVars(varargin{1}, mesh, vars, n_samples);
%     Q_app = getStateVars(varargin{2}, mesh, vars, n_samples);
%    % Get errors of the first kind (observation error)
%     varargout{1} = obj.get_errors(varargin{1}, varargin{2}, true);
%     % Get errors of the second kind (variable error)
%     varargout{2} = obj.get_errors(Q_ref', Q_app', true);
else
    % Difference matrix
    Delta = getError(varargin{1}, varargin{2}, 2, 0);
    % Scaling factors
    [~, ~, sf] = center_scale_data(varargin{1}, sf_C);
    % Scale the difference matrix
    Delta = Delta .* repmat(sf, 1, size(varargin{1},2)); % Multiply
    Delta = Delta ./ (repmat(sf.^2, 1, size(varargin{1},2)) + a_tol); % Divide
    % Rearrange difference matrix
    n_samples = size(varargin{1},2) / size(mesh,1);
    Y_delta = leaveStateVars(Delta, mesh, vars, n_samples);
    % First kind errors
    varargout{1} = mean(Y_delta, 1);
    % Second kind errors
    varargout{2} = mean(Delta, 2);
%     % Rearrange matrices
%     n_samples = size(varargin{1},2) / size(mesh,1);
%     Y_ref = leaveStateVars(varargin{1}, mesh, vars, n_samples);
%     Y_app = leaveStateVars(varargin{2}, mesh, vars, n_samples);
%     % Get errors of the first kind (observation error)
%     varargout{1} = obj.get_errors(Y_ref, Y_app);
%     % Get errors of the second kind (variable error)
%     varargout{2} = obj.get_errors(varargin{1}', varargin{2}', true);
end

end

