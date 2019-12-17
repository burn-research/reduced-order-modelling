function [X_noised] = add_noise_to_data(X, snr, remove_negative_terms)
% This function adds random noise to data.
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

    noised_column = awgn(X(:,col), snr, 'measured');

    % Remove negative terms if needed:
    if remove_negative_terms == true
        noised_column(noised_column<0) = 0;
    end

    X_noised = [X_noised, noised_column];

end
