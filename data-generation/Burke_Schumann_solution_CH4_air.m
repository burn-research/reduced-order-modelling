function [Z, state_space, state_space_labels, s] = Burke_Schumann_solution_CH4_air(n_points, graph)
% This function generates the thermo-chemical state-space from the
% Burke-Schumann solution. The fuel is pure CH4 and the oxidizer stream is
% air. The function returns the state_space matrix of size (n_points x 6)
% where the columns represent respectively: T, CH4, O2, CO2, H2O, N2.
%
% Input:
% ------------
% - n_points
%         number of points for mixture fraction sampling.
%
% - graph
%         boolean specifying whether you want to plot and save a graph of
%         the thermo-chemical state-space.
%
% Output:
% ------------
% - Z
%         mixture fraction vector.
%
% - state_space
%         matrix of the thermo-chemical state-space.
%
%              T  CH4  O2  CO2  H2O  N2
%             [                        ]
%             [                        ]
%             [                        ]
%             [                        ]
%             [                        ]
%             [                        ]
%             [                        ]
%
% - state_space_labels
%         string labels for every column.
%
% - s
%         stoichiometric ratio.

Tf_st = 2263.3; % stoichimetric flame temperature for CH4

% Mass fractions in (1) fuel stream and (2) oxidizer stream:
yCH4_1 = 1; % pure CH4
yO2_2 = 0.232; % air
yN2_2 = 1-yO2_2; % air

% Number of carbon atoms nC and hydrogen atoms nH in methane:
nC = 1;
nH = 4;

% Molecular weights:
MW_C = 12.0107/1000; % kg/mole
MW_H = 1.00784/1000; % kg/mole
MW_O = 15.999/1000; % kg/mole
MW_CH4 = MW_C + 4*MW_H; % kg/mole
MW_CO2 = MW_C + 2*MW_O; % kg/mole
MW_H2O = 2*MW_H + MW_O; % kg/mole
MW_O2 = 2*MW_O; % kg/mole

% Thermodynamic properties:
T1 = 300; %K
T2 = 300; %K

Q = 50.0*10^6*MW_CH4; % heat of combustion (LHV) [J/mole]
Cp = 1.4*1000; % heat capacity [J/kg-K]

% Stoichiometric ratio:
s = 2*MW_O2/MW_CH4;

% Stoichiometric mixture fraction:
Z_st = yO2_2 / (s*yCH4_1 + yO2_2)

% Mixture fraction vector:
Z = linspace(0, 1, n_points)';

% Variable vectors:
T = zeros(size(Z));
yCH4 = zeros(size(Z));
yO2 = zeros(size(Z));
yCO2 = zeros(size(Z));
yH2O = zeros(size(Z));
yN2 = zeros(size(Z));

% CH4 mass fraction:
for i = 1:1:length(Z)
    if Z(i) >= Z_st
        yCH4(i) = yCH4_1 * ((Z(i) - Z_st)/(1 - Z_st));
    else
        yCH4(i) = 0;
    end
end

% O2 mass fraction:
for i = 1:1:length(Z)
    if Z(i) <= Z_st
        yO2(i) = yO2_2 * (1 - Z(i)/Z_st);
    else
        yO2(i) = 0;
    end
end

% CO2 mass fraction:
yCO2_st = yCH4_1 * Z_st * (nC * MW_CO2/MW_CH4);

for i = 1:1:length(Z)

    if Z(i) <= Z_st
        yCO2(i) = yCO2_st * Z(i)/Z_st;
    else
        yCO2(i) = yCO2_st * ((1 - Z(i))/(1 - Z_st));
    end

end

% H2O mass fraction:
yH2O_st = yCH4_1 * Z_st * (nH * MW_H2O/(2*MW_CH4));

for i = 1:1:length(Z)

    if Z(i) <= Z_st
        yH2O(i) = yH2O_st * Z(i)/Z_st;
    else
        yH2O(i) = yH2O_st * ((1 - Z(i))/(1 - Z_st));
    end

end

% N2 mass fraction:
for i = 1:1:length(Z)
    yN2(i) = 1 - yCH4(i) - yCO2(i) - yH2O(i) - yO2(i);
end

% Temperature:
for i = 1:1:length(Z)

    Tu = T2 + Z(i) * (T1 - T2);

    if Z(i) <= Z_st
        T(i) = Tu + Q*yCH4_1/(Cp * 1 * MW_CH4) * Z(i);
    else
        T(i) = Tu + Q*yO2_2/(Cp * 2 * MW_O2) * (1 - Z(i));
    end

end

state_space = [T, yCH4, yO2, yCO2, yH2O, yN2];
state_space_labels = {'T', 'CH4', 'O2', 'CO2', 'H2O', 'N2'};

%% Plot state-space:
if graph == true

    figure();
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.7, 1]);
    subplot(2, 3, 1)
    plot(Z, T, 'k', 'LineWidth', 2)
    % xlabel('F [-]', 'FontSize', 25)
    ylabel('T [K]', 'FontSize', 25)
    set(gca, 'FontSize', 20); box on; grid on
    subplot(2, 3, 2)
    plot(Z, yCH4, 'k', 'LineWidth', 2)
    % xlabel('F [-]', 'FontSize', 25)
    ylabel('$Y_{CH4}$ [-]', 'FontSize', 25)
    set(gca, 'FontSize', 20); box on; grid on
    subplot(2, 3, 3)
    plot(Z, yO2, 'k', 'LineWidth', 2)
    % xlabel('F [-]', 'FontSize', 25)
    ylabel('$Y_{O2}$ [-]', 'FontSize', 25)
    set(gca, 'FontSize', 20); box on; grid on
    subplot(2, 3, 4)
    plot(Z, yCO2, 'k', 'LineWidth', 2)
    % xlabel('F [-]', 'FontSize', 25)
    ylabel('$Y_{CO2}$ [-]', 'FontSize', 25)
    set(gca, 'FontSize', 20); box on; grid on
    subplot(2, 3, 5)
    plot(Z, yH2O, 'k', 'LineWidth', 2)
    % xlabel('F [-]', 'FontSize', 25)
    ylabel('$Y_{H2O}$ [-]', 'FontSize', 25)
    set(gca, 'FontSize', 20); box on; grid on
    subplot(2, 3, 6)
    plot(Z, yN2, 'k', 'LineWidth', 2)
    % xlabel('F [-]', 'FontSize', 25)
    ylabel('$Y_{N2}$ [-]', 'FontSize', 25)
    set(gca, 'FontSize', 20); box on; grid on

    filename = ['Burke_Schumann_CH4_air_state_space.png'];
    saveas(gcf, filename, 'png');

end

end
