function plot_cluster_populations(obj, varargin)

pop = obj.getClusterPopulations();
bar(1:length(pop), pop); grid on;
title('Cluster populations')

end


