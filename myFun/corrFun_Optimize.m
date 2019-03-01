function [I, nice_guess] = corrFun_Optimize(samples, values, prediction_points, trendFun, corrFun_list, varargin)
% Input: if test points are present, they will be used to choose the best
% correlation function
test_data = [];
if ~isempty(varargin)
    test_data = varargin{1}; test_data = []; % Set off
end

% Initialize
mse_check = zeros( length(corrFun_list), 1 ); % MSE returned by Kriging
guesses = [0, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1]; % List of guesses to try
mse_guess = zeros( length(guesses), 1 ); % Means of the MSE returned for that guess
mse_saved = cell( length(guesses), 1 ); % Complete MSE returned for that guess
j_saved = zeros( length(guesses), 1 ); % Indeces of the guess that worked best for each corrFun

% The best guess is the one that returned the smallest mse
for i = 1 : length(corrFun_list) 
  % Run it several times in order to find the best GUESS for the hyperparameters
  val_kriged = cell(length(guesses), 1);
    parfor j = 1 : length(guesses)
        try
            % Perform Kriging
            [val_kriged{j}, mse_saved{j}] = performKriging(...
                samples, values, prediction_points,...
                trendFun, corrFun_list{i}, guesses(j)...
                );
        catch ME
            % In case the algorithm fails, return bad outputs
            val_kriged{j} = values(:,1) * 0 + Inf; 
            mse_saved{j} = values(:,1) * 0 + Inf;
        end
        
        % Save info about how things went for each guess
        if isempty(test_data)
            % MSE check
            temp = abs(mse_saved{j}) ./ (abs(val_kriged{j}) + eps);
            mse_guess(j) = metric_for_mse(temp);
        else
            % Test points check
            temp = abs(val_kriged{j} - test_data) ./ (abs(test_data) + eps);
            mse_guess(j) = metric_for_test(temp);
        end
    end
    
    % Index of the starting guess that worked best
    [~, j] = min(mse_guess);    
    mse = mse_saved{j}; % Get the corresponding MSE
    vals = val_kriged{j}; % Get the corresponding predictions

    % Store the best starting guess for each corrFun that's run
    j_saved(i) = j; 
    
    % Save info about how things went for each corrFun
    if isempty(test_data)
        % MSE check
        temp = abs(mse) ./ (abs(val_kriged{j}) + eps);
        mse_check(i) = metric_for_mse(temp);
    else
        % Test points check
        temp = abs(vals - test_data) ./ (abs(test_data) + eps);
        mse_check(i) = metric_for_test(temp);
    end
end

% Index of the corrFun that worked best
[~, I] = min(mse_check);

% Guess that worked best for that corrFun
nice_guess = guesses(j_saved(I)); % The corresponding best guess

end



function y = metric_for_mse(temp)

y = max( mean( temp, 1) );

end

function y = metric_for_test(temp)

y = max( mean( temp, 1) );

end


