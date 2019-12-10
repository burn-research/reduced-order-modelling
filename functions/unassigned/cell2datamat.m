function Y = cell2datamat(A)

% Number of cells
l = length(A);

Y = [];

for i = 1 : l
    Y = [Y, A{i}'];
end



end

