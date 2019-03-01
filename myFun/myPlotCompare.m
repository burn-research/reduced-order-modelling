function [varargout] = myPlotCompare(varargin)
%% Description
% Tired of habing to type something like:
%   plot(x1); hold on; plot(x2); hold off; grid on;
% each time I had to compare two vectors, here I write a small function
% that does so.
% ------------------------------------------------------------------------

%% Input processing

% x-vector: vector for the x-axis;
% y-vector: vector for the y-axis.

% User must provide 2 y-vectors at least
if nargin < 2
    error('Provide at least two vectors.');
end

% User is able to add x-vectors by doing "myPlotCompare(..., 'x', true)"
input_rule = find_string(varargin, 'x');
y_axis_log = find_string(varargin, 'y-axis');

% If y_axis_log is true, there are 2 input not to be plotted
if y_axis_log
    n_plots = nargin - 2;
else
    n_plots = nargin;
end

% If x-vectors are passed in, take off the input "x-axis" and half the number of plots
if input_rule
    n_plots = (n_plots - 2) / 2;
end

% Return an error if the number of provided vectors is not even (property
% 'x' is ON).
if input_rule && mod(nargin, 2) ~= 0
    error('Provide an even number of inputs when field "x" is ON. ');
end


%% Main
ii = 1;
while ii <= n_plots
    if ~input_rule
        if y_axis_log
            semilogy(varargin{ii})
        else
            plot(varargin{ii});
        end
        hold on;
        ii = ii + 1;
    else
        if y_axis_log
            semilogy(varargin{ii}, varargin{ii + 1});
        else
            plot(varargin{ii}, varargin{ii + 1});
        end
        hold on;
        ii = ii + 2;
    end
end
grid on; hold off;


end


function out = find_string(list_of_inputs, input_string)

out = false;

n_inputs = length(list_of_inputs);
for ii = 1 : n_inputs
    temp = list_of_inputs{ii};
    if ischar(temp) && strcmp(temp, input_string)
        out = list_of_inputs{ii+1};
    end
end

if ~islogical(out)
    error(['The value for the field "', input_string,'" must be logical.']);
end

end


