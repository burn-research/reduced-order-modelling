function [varargout] = xml2map(varargin)
%% DESCRIPTION


%% Inputs
nin = length(varargin); % Number of inputs

if nin < 1
    error('Not enough input arguments.');
end

% Input 1
keys = varargin{1};
if ~isa(keys,'cell')
    keys = num2cell(keys);
end

% Input 2
filename = ['/Users/gianmarcoaversano/opensmoke++suite-0.3.macos/',...
    'premixedLaminarFlame1D/ER_FC_T/T_298/']; % Default value
if nin > 1
    if ischar(varargin{2})
        filename = varargin{2};
    end
end


%% Main
data = xml2ws_levelup(filename,[]);
data_map = containers.Map(keys, data);


%% Output
varargout{1} = data_map;


end
