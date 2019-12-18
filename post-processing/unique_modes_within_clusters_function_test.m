clc, clear, close all

k = 2;
alpha = 0.99;

c1 = [1;2;1;2]/norm([1;2;1;2]);
c2 = [2;3;3;2]/norm([2;3;3;2]);
c3 = [3;5;3;1]/norm([3;5;3;1]);
c4 = [0;4;2;4]/norm([0;4;2;4]);
c5 = [5;5;5;5]/norm([5;5;5;5]);
c6 = [3;4;6;6]/norm([3;4;6;6]);
c7 = [7;0;7;7]/norm([7;0;7;7]);

modes_clusters = [c1, c6, c2, c3, c4, c2, c3, c4, c5, c5, c7, c1];

[unique_modes, no_unique] = unique_modes_within_clusters(modes_clusters, k, alpha);

clearvars -except modes_clusters unique_modes no_unique
