% Input: x - vector of integers
% Output: y - vector of integers without holes
function y=renumber(xin)

x=xin;
c=1:max(x);

for i=1:length(c)
    while ~ismember(c(i),x) && c(i)<max(x)
        x(x>c(i)) = x(x>c(i)) - 1;
    end
end

y=x;
end


