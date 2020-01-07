function [X_noised] = add_adaptive_noise_to_data(X, snr, remove_negative_terms)
% This function adds random noise to data. The noise is adapted to the
% absolute
%
% Input:
% ------------
% - X
%         the raw data set.
%
% - snr
%         Signal-to-Noise ratio. This parameter can be tuned to adjust the level of random noise.
%
% - remove_negative_terms
%         boolean specifying if negative terms in the new raw data matrix should be removed.
%         Possible use case: when the added noise make certain values negative but it is not physically desired.
%
% Output:
% ------------
% - X_noised
%         the raw data set with added noise.

%% add_noise_to_data()
[n_obs, n_vars] = size(X);

X_noised = [];

for col = 1:1:n_vars

    max_col = abs(max(X(:,col)));
    noised_column = zeros(size(X(:,col)));

    for i = 1:1:n_obs
        ratio = (2 - X(i,col)/max_col);
        noised_column(i) = awgn(X(i,col), ratio*snr, 'measured');
    end

    % Remove negative terms if needed:
    if remove_negative_terms == true
        noised_column(noised_column<0) = 0;
    end

    X_noised = [X_noised, noised_column];

end
