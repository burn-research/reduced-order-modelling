classdef DoE < handle
%% Description
%{
After creating the function fakeDoE() and it subfunctions, I've decided to
create a simple class that might help use fakeDoE().
%}


%% Properties
% Private
properties (Access = public)
    p = [];
    Yas = [];
    p_miss = [];
    Pot = [];
    Infl_rel = [];
    a = [];             % Number of samples / Total number of observations
end

% Public
properties (Access = public)
    p_orig = [];
    Y_orig = [];
    p0 = [];
    nSamplesMax = [];
    lhs = false;
end

% Hidden
properties (Hidden)
    cdp_filter = false;
end



%% Methods
% Constructor
methods
    function obj = DoE(varargin)
        
    end
end


% Set methods
methods
    function set.a(obj, val)
        if isempty(val)
            obj.a = [];
        elseif val >= 0 && val <= 1
            obj.a = val;
        else
            obj.a = .65;
            fprintf('\nWrong input provided.\n');
        end
    end
end


% Object methods
methods
    runDoE(obj); % Run fakeDoE() on this class object's properties
    getNumberOfSamples(obj, varargin); % Evaluates nSamplesMax 
    chooseP0(obj, n); % Chooses the samples that are the farthest as initial samples
    plotSamples(obj); % Plots the sampled input space
    plot_Influence(obj, varargin); % Plots the sampled input space with influences
    clear(obj); % Clears 
    [varargout] = createKpcaObj(obj, varargin);
end

    
end

