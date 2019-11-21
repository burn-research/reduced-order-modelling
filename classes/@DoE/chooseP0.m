function chooseP0(obj, n)

if n >= obj.nSamplesMax
    error('Provided input is too high.');
end

p = obj.p_orig;

% Center the points
[z, mu, ~] = zscore(p, 0, 2);
[~, l] = size(p);

% Evaluate distances from the mean
norm_of_z = sum(z.^2, 2).^.5;

% Sort distances
[z_sorted, I] = sort(norm_of_z, 'descend');

% Set p0
obj.p0 = p(I(1:n), :);
obj.p0 = sortrows(obj.p0);

end