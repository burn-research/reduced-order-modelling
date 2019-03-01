function [xdata, ydata, zdata] = getDataFromFigure(h, varargin)
% Get a handle to the current figure:
% h = gcf; %current figure handle
% The data that is plotted is usually a 'child' of the Axes object. 
% The axes objects are themselves children of the figure. 
% You can go down their hierarchy as follows:
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% Extract values from the dataObjs of your choice. 
% You can check their type by typing:
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
% Lines of code similar to the following would be required to bring the 
% data to MATLAB Workspace:
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
zdata = get(dataObjs, 'ZData');
end
