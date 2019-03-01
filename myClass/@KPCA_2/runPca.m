function runPca(obj, varargin)
%% Inputs
n_args = length(varargin);
if n_args > 0
    try
        obj.pca_approximation_order = varargin{1};
    catch ME
        disp(ME);
        return;
    end
end

%% Run PCA
if obj.is_local && ~obj.is_mesh_variable
    obj.centered_scaled_data = obj.subtractCentroid();
end

% Run PCA
obj.runPca_hiddenFunction();

% Useful
try
    obj.get_Kriging_targets();
catch 
    fprintf('- Could not get_Kriging_targets().\n');
end

% Recover data
obj.recoverPca();

% Estimate errors
try
    if ~obj.is_local
        obj.getPcaErrors();
    end
catch ME
    fprintf('[runPca:] Could not getPcaErrors()');
    disp(ME)
    for ii = 1 : length(ME.stack)
        disp(ME.stack(ii));
    end
end

end





