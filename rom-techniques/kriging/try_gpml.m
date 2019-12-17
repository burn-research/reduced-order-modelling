function [mu, s2, hyp2] = try_gpml(x, y, xs, varargin)

% Get useful sizes
y_dim = size(y,2);
n_points = size(xs,1);
x_dim = size(x,2);

% If more than one scalar output, use recursive function
if y_dim > 1
    % Initialize function outputs
    mu = zeros(n_points, y_dim);
    s2 = zeros(n_points, y_dim);
    hyp2 = cell(y_dim,1);
    % Rerun this same function, one scalar per run
    parfor ii = 1 : y_dim
        [mu(:,ii), s2(:,ii), hyp2{ii}] = try_gpml(x, y(:,ii), xs);
    end
    % Do not forget to stop here
    return
end

% Setup important regression properties
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; % Leave empty for no mean function
covfunc = @covSEiso; % Covariance function
likfunc = @likGauss; % Likelihood
infmethod = @infGaussLik; % Inference method

% Initialize the hyperparameters
hyp = struct('mean', [ones(x_dim,1); 1], ...
    'cov', [0 0], ...
    'lik', -1);

% Get best set of hyperparameters
hyp2 = minimize(hyp, @gp, -100, infmethod, meanfunc, covfunc, likfunc, x, y);

% Run the GP
[mu, s2] = gp(hyp2, infmethod, meanfunc, covfunc, likfunc, x, y, xs);

% Plot
% if x_dim == 1
%     figure()
%     f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
%     fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
%     hold on; plot(xs, mu); plot(x, y, '+'); hold off;
% end

end








