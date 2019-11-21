function y = aveVar(Y, vars, mesh)

[n, m] = size(Y);
n_var = length(vars);
n_points = size(mesh,1);

if (n_var ~= n/n_points)
    y = Y; return
end

y = size(n_var, m);
for i = 1 : n_var
    y_var = Y(1+n_points*(i-1):n_points*i,:);
    y(i,:) = mean( y_var , 1 );
end

end 

