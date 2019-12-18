% Artificial modes for a DNS1 data set
% The weights on each modes are set to 1.0 for a particular group of species:
%
%   [H2, O2, O, OH, H2O, H, HO2, H2O2, CO, CO2, HCO, N2]

% Assign weights on modes:
group_fuel = [1,0,0,0,0,0,0,0,1,0,0,0];
group_oxi = [0,1,0,0,0,0,0,0,0,0,0,1];
group_light = [0,0,1,1,0,0,0,0,0,0,0,0];
group_h = [0,0,0,0,0,1,0,0,0,0,0,0];
group_prod = [0,0,0,0,1,0,0,0,0,1,0,0];
group_hco = [0,0,0,0,0,0,0,0,0,0,1,0];
group_hxox= [0,0,0,0,0,0,1,1,0,0,0,0];

% Normalize the modes:
fuel = group_fuel/norm(group_fuel);
light = group_light/norm(group_light);
h = group_h/norm(group_h);
prod = group_prod/norm(group_prod);
hco = group_hco/norm(group_hco);
hxox = group_hxox/norm(group_hxox);
oxi = group_oxi/norm(group_oxi);

% Combine into modes matrix:
modes = [fuel' oxi' light' h' hxox' hco' prod'];
