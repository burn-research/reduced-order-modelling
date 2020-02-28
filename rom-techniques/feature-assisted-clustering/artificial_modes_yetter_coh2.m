% Artificial modes for chemical mechanism for burning CO/H2 [1].
%
% The mechanism contains 12 species:
%
% [H2, O2, O, OH, H2O, H, HO2, H2O2, CO, CO2, HCO, N2]
%
% The modes are constructed such that a particular species or a group of species
% are either "on" or "off". Modes are then normalized.
%
% References:
% -----------
% [1] Yetter, R. A., Dryer, F. L., & Rabitz, H. (1991). A comprehensive reaction
% mechanism for carbon monoxide/hydrogen/oxygen kinetics. Combustion Science and
% Technology, 79(1-3), 97-128.
%
% [2] Sutherland, J. C., Smith, P. J., & Chen, J. H. (2007). A quantitative
% method for a priori evaluation of combustion reaction models. Combustion
% Theory and Modelling, 11(2), 287-303.

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
