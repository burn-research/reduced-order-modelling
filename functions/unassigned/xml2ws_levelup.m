function [ varargout ] = xml2ws_levelup( inp_filename, second_level_name, varargin )
%XML2WS_LEVELUP Summary of this function goes here
%   Detailed explanation goes here


%% Default variables
default_inp_filename = ['/Users/gianmarcoaversano/opensmoke++suite-0.3.macos/',...
    'premixedLaminarFlame1D/ER_FC_T/T_298/'];
default_second_level_name = {'Output_04/', 'Output_05/', 'Output_06/',... 
    'Output_07/', 'Output_08/', 'Output_09/', 'Output_10/'};
numCases = 21; 


%% Inputs

% MANDATORY

% inp_filename
if (ischar(inp_filename) && strcmp(inp_filename, 'default')) || isempty(inp_filename)
    inp_filename = default_inp_filename;
elseif ~ischar(inp_filename) && ~isempty(inp_filename)
    error('1st input must be of type CHAR or empty.');
end

% second_level_name
if ischar(second_level_name) && strcmp(second_level_name, 'default') || isempty(second_level_name)
    second_level_name = default_second_level_name;
elseif ~ischar(second_level_name) && ~isempty(second_level_name)
    error('2nd input must be of type CHAR or empty.');
end


% OPTIONAL
nin = length(varargin); % Number of optional inputs provided

% 1-st input: numCases
if nin > 0
    numCases = varargin{1};
end


%% Main
l = length(second_level_name); % Number of folders
data = cell(l, 1);             % Initialize DATA

% Use xml2workspace(), output is stored in the cell array DATA
for i = 1 : l
    filename = [inp_filename, second_level_name{i}]; % Create FILENAME
    data{i} = xml2workspace(filename, numCases);
end


%% Output
varargout{1} = data;


end





