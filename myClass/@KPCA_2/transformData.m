function [varargout] = transformData(obj, fun, varargin)
%% Description:
% Training and original data are transformed (e.g. log-transformed).
%
% INPUT:
% - fun: (scalar) function to be used for the transformation.
%

%% Input
% Check input FUN is a function handle
isfun = @(f) isa(f, 'function_handle'); % Create your version for isfun()
if ~isfun(fun)
    error('Please provide a FUNCTION_HANDLE as input.');
end
% Input info
in_count = 0;
n_args = length(varargin);
% First varargin: if TRUE, predictions will also be transformed
if n_args > in_count
    in_count = in_count + 1;
    transform_predictions = varargin{in_count};
else
    transform_predictions = false;
end

%% Main
% Transform training data
obj.training_data = fun(obj.training_data);
% Transform original data
obj.original_data = fun(obj.original_data);
% Transform predictions
if transform_predictions
    obj.pca_recovered_data = fun(obj.pca_recovered_data);
    obj.cpca_recovered_data = fun(obj.cpca_recovered_data);
    obj.local_recovered_data_pca = fun(obj.local_recovered_data_pca);
    obj.local_recovered_data_cpca = fun(obj.local_recovered_data_cpca);
    obj.kriged_direct_data = fun(obj.kriged_direct_data);
    obj.kpca_predictions = fun(obj.kpca_predictions);
    obj.kcpca_predictions = fun(obj.kcpca_predictions);
    obj.klpca_predictions = fun(obj.klpca_predictions);
    obj.klcpca_predictions = fun(obj.klcpca_predictions);
end

%% Output
out_count = 0;
if nargout > out_count
    out_count = out_count + 1;
    varargout{out_count} = [];
end

end


