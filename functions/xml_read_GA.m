function [varargout] = xml_read_GA(filename)

xmlfile = fullfile(filename);
    
Pref.Str2Num = 'always';    
[tree, RootName, DOMnode] = xml_read(xmlfile, Pref); % Downloaded toolbox

% tree.profiles: matrix, columns are one var profile along the x-axis

varargout{1} = tree.profiles; % Saving data into matrices
varargout{2} = tree;

end