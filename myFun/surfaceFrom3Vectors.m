function surfaceFrom3Vectors(x, y, z)

f = fit([x,y],z, 'lowess');
plot(f);

end

