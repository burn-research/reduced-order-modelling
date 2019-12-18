clc, clear, close all

k = 3;
alpha = 0.95;

g1 = [1;2;1;2]/norm([1;2;1;2]);
g2 = [2;3;3;2]/norm([2;3;3;2]);
g3 = [3;5;3;1]/norm([3;5;3;1]);
g4 = [0;4;2;4]/norm([0;4;2;4]);

l1 = [5;5;5;5]/norm([5;5;5;5]);
l2 = [3;4;6;6]/norm([3;4;6;6]);
l3 = [7;0;7;7]/norm([7;0;7;7]);

modes_global = [g1, g2, g3, g4, g2, g3, g4];
modes_local = [l1, l2, l3, g2, g2, l1];

[rep_modes_global, rep_modes_local, cluster_annotation, no_repeating] = repeating_modes_global(modes_global, modes_local, k, alpha);

clearvars -except modes_global modes_local repeating_modes_global repeating_modes_local cluster_annotation no_repeating
