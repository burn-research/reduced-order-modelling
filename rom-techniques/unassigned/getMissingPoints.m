function p_miss = getMissingPoints(inp1, inp2)
%% Description



%% Input
% Chech which is the longest
if size(inp1,1) > size(inp2,1)
    p1 = inp1;
    p2 = inp2;
else
    p1 = inp2;
    p2 = inp1;
end

% Sizes
[r1, c1] = size(p1);
[r2, c2] = size(p2);

% Check size consistency
if c1 ~= c2
    error('Vectors must have same number of columns.');
end


%% Main
% Initialize logic vector
I = false(r1,1);

% Check if p1(i) is in p2
for i = 1 : r1
    for j = 1 : r2
        temp = prod( p1(i,:) == p2(j,:) ) == 1;
        if temp
            I(i) = true;
        end
    end
end


%% Output
p_miss = p1(~I,:);


end


