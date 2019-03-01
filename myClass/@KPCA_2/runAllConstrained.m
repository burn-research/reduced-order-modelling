function runAllConstrained(obj, varargin)

% Run all routines
obj.runCpca();
obj.runCpcaKriging();
obj.runClustering();
obj.runLcpca();
obj.runLcpcaKriging();

end


