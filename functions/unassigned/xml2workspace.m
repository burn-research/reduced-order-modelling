function [varargout] = xml2workspace(inp_filename, varargin)
%% DESCRIPTION
%{

%}

%% Default variables
default_filename = ['/Users/gianmarcoaversano/opensmoke++suite-0.3.macos/',...
    'premixedLaminarFlame1D/ER_FC_T/T_298/Output_10/'];
numCases = 21; 
isPlot = false;
j1 = 1;
var1 = 'temperature';


%% Inputs
nin = length(varargin);
if nin > 0
    numCases = varargin{1};
end
if nin > 1
    isPlot = varargin{2}; % Plot: True / False
end
if nin > 2
    j1 = varargin{3}; % Which case to plot
end
if nin > 3
    var1 = varargin{4};
    if ~isa(var1, 'char')
        var1 = 'temperature';
    end
end


%% Main

% Getting the data from the Output.xml files
data = cell(numCases,1);
for i = 1 : numCases
    temp2 = ['Case',num2str(i-1),'/Output.xml'];
    filename = [inp_filename, temp2];
    
    [data_temp{i}, tree] = xml_read_GA(filename); % Saving data into matrices
end


temp = strsplit(tree.mass_DASH_fractions)';
indexSpecies = temp(2:3:end);
temp = strsplit(tree.additional)';
indexQuantities = [temp(2:3:end), temp(3:3:end)]; 
indexAll = [temp(2:3:end); indexSpecies];


t = 1:numCases; t = t(:); 
temp = cell(length(t), 1);
for i = 1 : length(t)
    temp{i} = num2str(t(i)); 
    temp{i} = strrep(temp{i}, '.', '');
end

for i = 1 : numCases
    data.case{i,1} = data_temp{i};
end

data.indexAll = indexAll; 
data.indexQuantities = indexQuantities; 
data.indexSpecies = indexSpecies;
data.cases = temp;



%% Output
varargout{1} = data;



%% Plot
if isPlot 
    % Find the index of the variable to plot
    i1 = find(strcmp(var1, indexAll));  
    
    % Plot
    x = data.case{j1}(:,1);
    y = data.case{j1}(:,i1);
    plot(x,y);
    xlabel('x [cm]'); 
    grid on; 
    legend(var1);
end


end



