function plotTrainingErrors(obj, varargin)

y2plot1 = 100 * obj.Err_PCA_obs;
y2plot2 = 100 * obj.Err_LPCA_obs;

clustNum = length(obj.lpca);

t = 'Training points';

dummyFun(y2plot1, y2plot2, clustNum, t);

end

function dummyFun(y2plot1, y2plot2, clustNum, t)
% INPUTS
%   y2plot1,2: values to plot
%   clustNum: number of clusters

setFigOpts;
figure();
semilogy(y2plot1, '--o'); grid on; hold on;
semilogy(y2plot2, '--o');

xlabel('Observation number [-]'); ylabel('Error [%]');
title(t); 
legend('K-PCA',['K-LPCA',clustNum]);


end