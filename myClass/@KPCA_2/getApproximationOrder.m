function k = getApproximationOrder(obj, treshold)

% Get approximation order
t = cumsum(obj.pca_eigenvalues) / sum(obj.pca_eigenvalues); % Work with [-]
i = find(t > treshold);
k = i(1);
this = num2str(round(100*t(k),2));
fprintf([this,' per cent of the data variance captured.\n\n']);

end
