function internalPlotFun(y2plot1, y2plot2, clustNum, t)

setFigOpts;
figure();
semilogy(y2plot1, '--o'); grid on; hold on;
semilogy(y2plot2, '--o');

xlabel('Observation number [-]'); ylabel('Error [%]');
title(t); 
legend('K-PCA',['K-LPCA',clustNum]);


end

