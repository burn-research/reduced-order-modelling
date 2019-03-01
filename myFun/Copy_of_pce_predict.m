function [ys] = pce_predict(H, c, xs, y)
if ~exist('xs', 'var') || isempty(xs)
    ys = [];
    return
end
if ~exist('y', 'var') || isempty(y)
    y = [];
end
% Determine 
if ~iscell(H)
    univariate = true;
else
    univariate = false;
end
% Predict
if univariate
    % Univariate
    ys = H(xs,:) * c;
else
    % Multi-variate
    x_dim = size(xs, 2);
    if length(H) ~= x_dim
        error('Length of H must equal x_dim.');
    end
    % Product of univariate polynomials 
    Hs = [];
    for i = 1 : x_dim - 1
        for j = i + 1 : x_dim
            for k = 1 : size(H{i}, 2)
                for l = 1 : size(H{j}, 2)
                    h = H{i}(xs(:,i),k) .* H{j}(xs(:,j),l);
                    Hs = [Hs, h];
                end
            end
        end
    end
    % Determination of coefficients
    ys = Hs * c;
end
% Post-processing
if ~isempty(y)
    [y0, mu] = center(y, 1);
    [~, sig] = scale(y, y0, 3);
    ys = unscale(ys, sig);
    ys = uncenter(ys, mu);
end
end

