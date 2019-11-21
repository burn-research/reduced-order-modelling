function setScores(obj, varargin)
%% Input
possible_models = {'pca', 'cpca', 'lpca', 'lcpca'};

n_args = length(varargin);

scores = [];
model = [];
kriged = false;
if n_args > 0
    % Get the scores 
    scores = varargin{1};
end
if n_args > 1
    % Get the model type
    model = varargin{2};
end
if n_args > 2
    % Flag for kriged scores
    kriged = varargin{3};
end


%% Main
switch model
    case possible_models{1} % PCA
        if kriged
            % Set the kriged scores
            obj.kriged_pca_scores = scores;
        else
            % Set the PCA scores
            obj.pca_scores = scores;
        end
    case possible_models{2} % CPCA 
        if kriged
            % Set the kriged scores
            obj.kriged_cpca_scores = scores;
        else
            % Set the CPCA scores
            obj.cpca_scores = scores;
        end
    case possible_models{3} % LPCA
        if ~iscell(scores)
            % A cell-array is needed in this case 
            error('Input (1) must be a cell array if model is LOCAL.');
        end
        if kriged
            % Set the kriged scores in every LPCA object
            for ii = 1 : length(scores)
                obj.local_pca{ii}.pca_scores = scores{ii};
            end
        else
            % Set the PCA scores in every LPCA object
            for ii = 1 : length(scores)
                obj.local_pca{ii}.kriged_pca_scores = scores{ii};
            end
        end
    case possible_models{4} % LCPCA
        if ~iscell(scores)
            % A cell-array is needed in this case
            error('Input (1) must be a cell array if model is LOCAL.');
        end
        if kriged
            % Set the kriged scores in every LPCA object
            for ii = 1 : length(scores)
                obj.local_pca{ii}.cpca_scores = scores{ii};
            end
        else
            % Set the CPCA scores in every LPCA object
            for ii = 1 : length(scores)
                obj.local_pca{ii}.kriged_cpca_scores = scores{ii};
            end
        end
    otherwise
        error('No possible model provided.')
end


end


