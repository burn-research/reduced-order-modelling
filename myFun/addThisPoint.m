function p_plus = addThisPoint(p_cdt, p0, I)
%% Description
% Adds one point from a set of points indicated by I, taken from p_orig, to
% p0. The first point in p_cdt is added if it is not already contained in
% p0, otherwise the second one is added, and so on.
%

%% Input
[ri, ci] = size(I);
[rp, cp] = size(p_cdt);
[r0, c0] = size(p0);

% Check size consistency
if ri ~= rp
    error('Wrong size provided.');
elseif ci ~= 1
    error('DIM2 should be one.');
elseif cp ~= c0
    error('DIM2 of these variables should be the same.');
end


%% Main
cond = true; % While-loop breaker
i = 1; % First index
while (cond)
    new_point = p_cdt(I(i),:); % Get the point
    % Check the point does not belong to p0
    J = ~ismember(new_point, p0, 'rows');
    if J
        cond = false;
    else
        i = i + 1;
    end
    if i == size(p_cdt,1)
        i;
    end 
end


%% Output
p_plus = new_point(1,:);


end

