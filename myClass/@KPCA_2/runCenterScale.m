function runCenterScale(obj, varargin)
% To avoid conflicts with the 2 properties center_scale and
% scaling_criterion, if either of those suggests no center-scaling, no
% center-scaling will occur then. The property center_scale is checked
% first and is thus given priority

if obj.center_scale
    [obj.centered_scaled_data, obj.mean_column, obj.scaling_factors] = ...
        center_scale_data(obj.training_data, obj.scaling_criterion, ...
        obj.mean_column, obj.scaling_factors);
else
    [obj.centered_scaled_data, obj.mean_column, obj.scaling_factors] = ...
        center_scale_data(obj.training_data, 0, ...
        obj.mean_column, obj.scaling_factors);
end

end


