function getRidOfNegatives( obj )

tol = 1e-22;


if ~isa(obj.Y, 'cell')
    obj.Y(obj.Y<tol) = 0;
    return;
end


l = length(obj.Y);

for i = 1 : l
    temp = obj.Y{i};
    temp(temp < tol) = 0;
    obj.Y{i} = temp;
end


end

