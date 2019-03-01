function [varargout] = plotRelInfl(p, Infl, varargin)
%% Description
%{
NOT READY

This function plots the Relative Influence of Adaptive Sampling. Usually, P
is given as a vector (observations x dimension), making it harder to use
the function surf(). As the plot3() function might not be good, this
function preprocess the data in order to use the function surf() at the
end.

Makes sense only if P's second dimension's size is 2!
%}


%% Input
[obs, dim] = size(p);
if obs ~= size(Infl, 1)
    error('Inputs have diffent sizes.');
end
if dim ~= 2
    error('First input must be a (any x 2) vector.');
end


% P


% Infl


% Varargin{1}


% Varargin{2}



%% Main
scatter3(p(:,1), p(:,2), Infl, Infl*0+100, Infl, 'filled');
grid on;
xlabel('Parameter 1');
ylabel('Parameter 2');
zlabel('Relative Influence [-]')
title('Relative Influence of the samples');


%% Output



end

