function setFigOpts(varargin)
%% Description
% My personal small function to set my favourite plotting settings
%

%%
set(0,'defaultaxesfontsize', 18);
set(0,'defaulttextfontsize', 18);
set(0,'DefaultLineLineWidth', 1.4);
set(0,'DefaultLineMarkerSize', 9);
set(0,'defaulttextfontweight', 'Bold');
set(0,'DefaultTextInterpreter', 'LaTex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

end