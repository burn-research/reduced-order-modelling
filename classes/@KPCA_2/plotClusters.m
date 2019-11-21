function plotClusters(obj, varargin)

% Default variables
marker = 'x';


% Input
nin = length(varargin);
if nin > 0
    marker = varargin{1};
    if ~isa(marker, 'char')
        marker = 'x';
    end
end


% Clusters
clusteredPoints = partitionVQ(obj.xp_rows, obj.idx, false);


% Plot
setFigOpts;
figure(); 
for i = 1 : length(clusteredPoints)
    plot(clusteredPoints{i}(:,1), clusteredPoints{i}(:,2), marker);
    hold on;
    leg{i} = num2str(i);
end
grid on;
legend( leg(:) );

xlabel('x [-]');
ylabel('Variables');
title('Clusters graphical representation');


end


% subplot(2,1,1), PLOT(income)
% subplot(2,1,2), PLOT(outgo)



