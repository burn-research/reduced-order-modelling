function setRegressionMethod(obj, method, varargin)

% Find which model the user wants
for ii = 1 : length(obj.reg_m)
    for jj = 1 : length(obj.reg_m{ii})
        if strcmp(obj.reg_m{ii}{jj}, method)
            this = ii; 
            obj.reg_method = method;
            break
        end
    end
end

% Set the regression flags
if this == 1 || this > 3
    % If Kriging
    obj.is_kriging = true;
    obj.is_gpml = false;
    obj.is_linear = false;
    obj.stacking = false;
elseif this == 2
    % If MARS
    obj.is_kriging = false;
    obj.is_gpml = false;
    obj.is_linear = true;
    obj.stacking = false;
elseif this == 3
    % If Gpml
    obj.is_kriging = false;
    obj.is_gpml = true;
    obj.is_linear = false;
    obj.stacking = false;
end

% Do the same for the local clusters
if ~obj.is_local && ~isempty(obj.local_pca)
    for ii = 1 : length(obj.local_pca)
        obj.local_pca{ii}.setRegressionMethod(method);
    end
end

end




