classdef NumericalData < handle
% Summary of this class:
%{
Useful to preprocess data before using the class KPCA.

Outputs of simulations are usually provided in cell array, where each cell
is an M-by-V matrix. M: observations, V: variables.

In order to make things easier, and set and get information from these data
faster, this class defines some useful methods that helps you do that.
%}
    

 
properties (SetAccess = private)
    % Constructed
    Y = {};
    
    % Others
    Q_mat = [];
    Y_mat = [];
    xp = [];
end

properties (Access = public)
    X = {};
    P = {};
    vars = {};
end
    
properties (Dependent)
    numberOfDatasets;
end





methods
    % Constructor
    function obj = NumericalData(Y, varargin)
        % Data-sets
        obj.Y = Y;
        
        % Additional inputs
        nin = length(varargin);
        if nin > 0
            obj.X = varargin{1};
        end
        if nin > 1
            obj.P = varargin{2};
        end
        if nin > 2
            obj.vars = varargin{3};
        end
        
        
    end
end


% Get methods
methods
    function val = get.numberOfDatasets(obj)
        val = length(obj.Y);
    end
end


% Set methods
methods
    
    % Mesh
    function set.X(obj, X)
        
        % If X is DOUBLE: common mesh case 
        if isa(obj.Y,'cell') && isa(X,'double')
            
            % Check if all datasets do have a common mesh
            for i = 1 : obj.numberOfDatasets - 1
                temp = (size(obj.Y{i}) == size(obj.Y{i+1}));
                if ~temp(1) && ~temp(2)
                    error('Something is wrong.');
                end
            end
            
            % Check if X has the same size of the dataset(s)
            temp = size(obj.Y{1}) == size(X);
            if ~temp(1) && ~temp(2)
                error('Invalid input.');
            end
        end
        
        
        % If X is CELL: different meshes
        if isa(obj.Y,'cell') && isa(X,'cell') 
            
            % Check number of meshes == number of data-sets
            if obj.numberOfDatasets ~= length(X)
                error('Different number of datasets');
            end
            
            % Check internal dimension are equal
            for i = 1 : obj.numberOfDatasets
                if size(obj.Y{i},1) ~= size(obj.X{i},1) && size(obj.Y{i},2) ~= size(obj.X{i},2)
                    error('Something is wrong.');
                end
            end
        end
        
        
        
        % Everything was OK
        obj.X = X;
    end
    
    
    % Parameters
    function set.P(obj, P)
        
        % Check number of samples == numberOfDatasets 
        if size(P,1) ~= obj.numberOfDatasets
            error('Too many or a few samples.');
        end
        
        % Everything was OK
        obj.P = P;
    end
    
    
    % Variables
    function set.vars(obj, vars)
        
        % Check number of samples == numberOfDatasets 
        temp = size(obj.Y{1});
        if length(vars) ~= temp(1) && length(vars) ~= temp(2)
            error('Too many or a few variables.');
        end
        
        % Everything was OK
        obj.vars = vars;
    end
    
end



% Object methods:
methods
    generateMat(obj);
    generateXP(obj);
    getRidOfNegatives(obj);
end



% Static methods: 
methods (Static)


end
    




end

