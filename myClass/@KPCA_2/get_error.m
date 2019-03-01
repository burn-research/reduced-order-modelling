function y = get_error(obj, x_original, x_predicted, varargin)

if nargin < 3
    y = obj.get_errors(x_original, x_predicted);
else
    y = obj.get_errors(x_original, x_predicted, varargin(:));
end

end