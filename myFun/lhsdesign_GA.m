function X = lhsdesign_GA(n_samples, dim, up_b, lo_b)
% Check sizes
if length(up_b) ~= length(lo_b) || dim ~= length(up_b)
    error('Wrong sizes. Size(up_b) = %i; Size(lo_up) = %i; dim = %i \n', length(up_b), length(lo_b), dim);
end
% Run LHS on a length-1 hypercube
X = lhsdesign(n_samples, dim, 'smooth', 'on');
% Rescale the hypercube
for ii = 1 : dim
    X(:,ii) = lo_b(ii) + X(:,ii) * (up_b(ii) - lo_b(ii));
end

end




