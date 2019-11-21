function [A] = index_set(M, p, q)
if q <= eps
    q = eps;
end
A = permn(0:p, M);
I = sum(abs(A).^q, 2).^(1/q) > p; % Elements to get rid of
A(I,:) = [];
A = A + 1;
end

