function [varargout] = samplingInfluence_Kriging(X, Y, trendFun, corrFun, hyp_guess, is_slow)
%% Description
% Input:
% - X: samples, obs x var
% - Y: values, obs x var
% - trendFun: order of the Kriging trend function (0, 1, 2, ecc)
% - corrFun: Kriging kernel, function handle
%

%% Input
N = size(X,1); % Number of samples
y_dim = size(Y,2); % Output dimension
% Check sizes
if N ~= size(Y,1)
    error('Number of rows must coincide.');
end
% Scale Y
[~, mu, sig] = center_scale_data(Y', 1);
mu = mu';   sig = sig';
% Switch from order to string
trendFun = ['regpoly', num2str(trendFun)];

%% Main
a_tol = 1e-16; % Tolerance
% Slow or Fast LOO?
Infl_rel = zeros(N,1); % Initialize output
Infl_rel_y = cell(N,1);
if is_slow
    % Slow LOO
    for ii = 1 : N
        % Leave-one-out
        I = true(N,1);  I(ii) = false;
        % Get the Kriging models
        [~, ~, k_models] = performKriging(X(I,:), Y(I,:)', [], trendFun, corrFun, hyp_guess);
        % Calculate error for each target
        Infl_rel_y{ii} = zeros(1, y_dim);
        for jj = 1 : y_dim
            Infl_rel_y{ii}(jj) = k_models{jj}.predict(X(ii,:)) - Y(ii,jj);
        end
        % Calculate error for this fold
        Infl_rel_y{ii} = Infl_rel_y{ii} ./ (sig + a_tol);
        clear k_models
    end
else
    % Fast LOO
    % Get the Kriging models now
    evalc('[~, ~, k_models] = performKriging(X, transpose(Y), [], trendFun, corrFun, hyp_guess);');
    for ii = 1 : N
        % Leave-one-out
        I = true(N,1);  I(ii) = false;
        Infl_rel_y{ii} = zeros(1, y_dim);
        for jj = 1 : y_dim
            % The previously found h-parameters are kept
            evalc('k_models{jj} = k_models{jj}.fit(X(I,:), Y(I,jj));'); 
            % Calculate error for this fold, y by y 
            Infl_rel_y{ii}(jj) = k_models{jj}.predict(X(ii,:)) - Y(ii,jj);
            % Clean-up
            k_models{jj} = k_models{jj}.cleanup();
        end
        % Calculate error for this fold 
        Infl_rel_y{ii} = Infl_rel_y{ii} ./ (sig + a_tol);
    end
end
% Get influence of each point
parfor ii = 1 : N
    Infl_rel(ii) = norm(Infl_rel_y{ii});
end

%% Output
if nargout > 0
    varargout{1} = Infl_rel;
end

end

function Y = transpose(X)
Y = X';
end


