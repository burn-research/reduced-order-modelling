classdef input_check < handle
%% Description


%% Properties
properties (Access = public)
    value = [];
end
properties (SetAccess = private)
    type;
    size;
end

%% Constructor
methods
    function obj = input_check(variable, varargin)
        obj.value = variable;
        obj.type = class(variable);
        obj.size = size(variable);
    end
end

%% Set methods
methods
    function set.value(obj, val)
        if isempty(obj.value)
            obj.value = val; return;
        end
        if ~strcmp(class(val), obj.type)
            error('You cannot change this variable`s type.');
        end
        if any(size(val) ~= obj.size)
            error('You cannot change this variable`s size.');
        end
        obj.value = val;
    end
end

%% Get methods
methods
    
end

%% Methods
methods
    
end


end

