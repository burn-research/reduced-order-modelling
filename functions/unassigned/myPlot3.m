function myPlot3(h, varargin)

% Additional input
nin = length(varargin);

% 1st additional input
s = 'o';
if nin > 0
    s = varargin{1};
    if ~ischar(s)
        fprintf('\n This input should be a CHAR. \n'); s = 'o';
    end
end


% Input
[r, c] = size(h);
if r ~= 3 && c ~= 3
    error('Wrong input.');
elseif c == 3
    x = h;
elseif r == 3
    x = h';
end


% Main
setFigOpts; % Loading favorite options

figure(); % Plot
plot3(x(:,1), x(:,2), x(:,3), s);
grid on;


end

