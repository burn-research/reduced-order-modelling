function [p, p_miss, varargout] = discrete_lhs(p_orig, p0, nSamples)
%% Description
%{
p_orig, p0: samples x variables
%}


%% Input
% Sorting
[r, c] = size(p0); % c: input space dimension; r: number of points

% At least 3 training points have to be provided
if r < 1
    error('Not enough starting points.');
end

% Membership of p0 in p_orig
if any(~ismember(p0, p_orig))
    error('Asked to sample a point for which there is no data.');
end

% Get number of samples to be added
nSamplesToAdd = nSamples - r;


%% Main

% Run the continuous LHS
cdp = lhsdesign(nSamplesToAdd, c, 'smooth','off'); % Candidate points LHS

% Evaluate the upper and lower bounds
B = zeros(2, c);
for i = 1 : c
    B(1, i) = max(p_orig(:,i));
    B(2, i) = min(p_orig(:,i));
end

% Re-scale the LH
for i = 1 : size(cdp, 2) 
    delta = B(1,i) - B(2,i);    
    offset = B(2,i);
    cdp(:,i) = cdp(:,i) .* delta + offset;
end

% Get the closest discrete points
p = p0;
for i = 1 : size(cdp, 1)
    d = p_orig - repmat(cdp(i,:), size(p_orig,1), 1); % Evaluate distances
    d = ( sum(d.^2, 2) ).^.5; % Evaluate norm of distance
    [~, I] = sort(d); % Sort distances in ascending order
    try
        p = [p; addThisPoint(p_orig, p0, I)]; % Add this point
    catch
        p = p;
    end
end

%% Output
% Chosen points
p = sortrows(unique(p,'rows'));

% Discarded points
I = ~ismember(p_orig, p, 'rows');
p_miss = p_orig(I,:);

if nargout
    varargout{1} = ~I;
end


end

