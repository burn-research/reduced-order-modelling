function generateMat(obj)

obj.Q_mat = cell2datamat(obj.Y);

if isa(obj.X,'double')
    obj.Y_mat = leaveStateVars(obj.Q_mat, obj.X, length(obj.vars), size(obj.P,1));
end
    
    
end

