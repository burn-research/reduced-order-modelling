function xvalScore = runLOO(samples, values, tf, cf, ls)
%% Description
% LOO validation for the Kriging class.
% cvpe() does not re-evaluates the length-scales, so here I write a
% functions that does.
% The input 'ls' has to be some optimized length-scales. In this function,
% they will be used as starting guess for each new SM.
%

%% Main
n_samples = size(samples,1);
xval = zeros(n_samples,1);
% LOO-cross validation
for ii = 1 : n_samples
    I = true(n_samples,1); I(ii) = false;
    % Fit cross-validated Kriging model
    temp = values(I,:)';
    evalc('[~, ~, ki] = performKriging(samples(I,:), temp, [], tf, cf, ls);');
    % Calculate error for this fold
    xval(ii,:) = mean((ki{1}.predict(samples(ii,:)) - values(ii,:)).^2);
    clear ki
end
% Final score is the mean of all fold errors
xvalScore = mean(xval);
end





