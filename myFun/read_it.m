function [data] = read_it(filename, numHeaders, numDataLines, colNum)
% Create format
fmt = [ repmat('%*s', 1, colNum-1), '%f%[^\n]'];
% Open file
fid = fopen(filename, 'rt');
% Read text
data = textscan(fid, fmt, numDataLines, 'HeaderLines', numHeaders);
% Close file
fclose(fid);
end



