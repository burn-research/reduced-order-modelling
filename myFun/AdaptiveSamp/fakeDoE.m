function [p, Yas, p_miss, Pot, Infl_rel, I_sorted, Infl_log, X_log, varargout] = fakeDoE(x_orig, Y_orig, x0, nSamplesMax, varargin)
%% Description
%{
This function performs a fake Design of Experiment.
Given an 'original' dataset Y_orig, of M+W observations, this function
chooses the best M osbervations, with M provided by the user.
The term 'fake' is due to the fact that, once Adaptive Sampling has chosen
the new point to sample, the observation is added to the dataset, but this
observation is not the output of any simulation: we already have
the data!
This can be useful when performing validation on predicting models: we have
a dataset, we select the training observations (with the Adaptive Sampling
technique) and use the remainder for validation of the predictions.

INPUTS
p_orig:
Y_orig:
p0:
nSamplesMax:

OUTPUTS
p:
Yas:
p_miss;
Infle_rel:
%}

%% Inputs
n_start = size(x0,1);
% Size
[rows, cols] = size(Y_orig);
% Sorting
[r, c] = size(x0); % c: input space dimension; r: number of points
if c == 1
    x_orig = sort(x_orig);
    x0 = sort(x0);
end
% At least 3 training points have to be provided
if r < 3
    error('Not enough starting points.');
end
% Final number of samples (check)
if nSamplesMax > cols
    error('Too many samples demanded. nSampleMax > number of cols of Y_orig.');
end
% Consistency: p_orig and Y_orig
if cols ~= length(x_orig)
    error('Number of original points not equal to number of columns of Y_orig.');
end
% Membership of p0 in p_orig
if any(~ismember(x0, x_orig))
    error('Asked to sample a point for which there is no data.');
end
% Optional inputs
nin = length(varargin);
% Filter for CDP
cdp_filter = false;
if nin > 0
    cdp_filter = varargin{1};
end

%% Fake DoE
% Candidate points
cdp = x_orig;             
i = ismember(cdp, x0, 'rows');     
cdp(i,:) = [];
% The initial data set
Yas = sampleit(x0, x_orig, Y_orig);
Infl_log = {};
X_log = {};
Pot = [];
% Enrich the data set
while size(x0, 1) < nSamplesMax && ~isempty(cdp)
    % Get a new point
    tic;
    [newpoint_loc, pot, Infl_log{end+1}] = AdaptiveSampling(x0, Yas, cdp);
    newpoint = cdp(newpoint_loc,:);
    time_as = toc; fprintf('\nAS: %.2f s. ', time_as);
    % Update the training points
    x0 = [x0; newpoint]; % x0 = sortrows(x0); 
    X_log{end+1} = x0;
    % Save Pot and Infl
    Pot = [Pot; max(pot)]; 
    % Update the candidate points
    i = ismember(cdp, x0, 'rows'); 
    cdp(i,:) = [];
    % Update the data set
    Yas = sampleit(x0, x_orig, Y_orig);
    % Faster version
    if cdp_filter
        n_max = floor(.2 * (size(x_orig,1) - size(x0,1)));
        cdp = discrete_lhs(x_orig, x0, n_max);
    end
    % Print status
    prog = 100 * (length(x0) / nSamplesMax);
    fprintf('\nProgress is at %.2f per cent.\n', prog);
end
% Influences
Infl_rel = samplingInfluence_PCA(x0, Yas, 0, 0);
% Training points chosen
[p, I] = sortrows(x0);
Infl_rel = Infl_rel(I);
I_sorted = I - n_start;
% Get the prediction points
I = ~ismember(x_orig, p, 'rows');
p_miss = x_orig(I,:);

%% Outputs
% p, Yas, p_miss, Pot, Infl_rel

%% Varargout

end




