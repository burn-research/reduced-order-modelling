function [ortho_modes] = orthogonalize_factors(modes)
% This function performs Gram-Schmidt orthogonalization of the basis
% vectors (factors) stored as columns in the matrix `modes`.
%
% Input:
% ------------
% - modes
%       trimmed, normalized factors.
%
% Output:
% ------------
% - ortho_modes
%       trimmed, orthonormalized factors.

%% orthogonalize_factors()
[m, n] = size(modes);
ortho_modes = zeros(m, n);
R = zeros(n, n);

for i = 1:size(modes, 2)

    col = modes(:,i);

    for j = 1:i-1

        R(j,i) = ortho_modes(:,j)' * modes(:,i);
        col = col - R(j,i) * ortho_modes(:,j);

    end

    R(i,i) = norm(col);
    ortho_modes(:,i) = col/R(i,i);

end
