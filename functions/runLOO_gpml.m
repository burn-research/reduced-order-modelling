function xvalScore = runLOO_gpml(samples, values, gpr_opts, gpr_hyp, h_optimize)
%% Description
% LOO validation for the gpml toolbox.
% With or without h-parameters re-evaluation.
% The input has to be some optimized length-scales. In this function,
% they will be used as starting guess for each new SM.
% INPUT:
% - samples: obs x var
% - values: obs x var
%

%% Main
n_samples = size(samples,1);
xval = zeros(n_samples,1);
% LOO-cross validation
for ii = 1 : n_samples
    I = true(n_samples,1); I(ii) = false;
    % Fit cross-validated Kriging model
    evalc('[mug, sig, hyp] = runGpml(samples(I,:), values(I,:), samples(~I,:), gpr_opts, gpr_hyp, h_optimize);');
    % Calculate error for this fold
    xval(ii,:) = mean((mug - values(~I,:)).^2);
end
% Final score is the mean of all fold errors
xvalScore = mean(xval);
end




