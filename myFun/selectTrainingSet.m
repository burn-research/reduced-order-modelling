function [samples, values, I] = selectTrainingSet(...
            is_mesh_variable, dim_mesh, mesh_size, samples, values, n_total, smp) 

% Size conformity check
n = size(samples,1) / mesh_size; % Probably the number of variables
if n ~= round(n) 
    error('[update_Kriging:selectTrainingSet] This number should be an integer!');
end
        
% The user has already selected the best mesh points
if ~isempty(smp) % SMP: selected mesh points
    % To be developed
    I = ismember(samples(1:mesh_size,1:dim_mesh), smp, 'rows');
    I = repmat(I, 1, n);
    values = values(I,:);
    samples = smp;
    return
end

% We need to select a subset of the training samples that involve all the
% parameters but only a subset of the mesh points

% Get the starting mesh points randomly
n0 = 10; % Number of starting points
Is = sort( max(round( rand(n0,1) * mesh_size ), 1) );

% n_total is for samples: we need to choose a subset of grid points
num_grid_points = round(n_total / n);

% Run the discrete LHS rountine, choose a subset of mesh points
temp = samples(1:mesh_size,1:dim_mesh);
[~, ~, I] = discrete_lhs(temp, temp(Is,:), num_grid_points);        
I = repmat(I, n, 1); % Repeat the vector

% Output
values = values(:,I);
samples = samples(I,:);

end



