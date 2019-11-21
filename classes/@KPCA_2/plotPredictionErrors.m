function plotPredictionErrors(obj, varargin)

I = ismember(obj.xp_sorted, obj.xp_kriged, 'rows');
y2plot1 = 100 * obj.Err_KPCA_obs(I);
y2plot2 = 100 * obj.Err_KLPCA_obs(I);

clustNum = length(obj.lpca);

t = 'Prediction points';

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
