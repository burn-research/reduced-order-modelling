function decoded_data = uncenterUnscale(obj, decoded_data, mean_column, scaling_factors, varargin)

% Relavant dimension
N = size(decoded_data, 1); 

% Diagonal matrix of scaling factors
D = spdiags(scaling_factors, 0, N, N);

% Unscale
decoded_data = D * decoded_data;

% Uncenter
M = repmat(mean_column, 1, size(decoded_data, 2));

% Decoded data
decoded_data = M + decoded_data;


end

