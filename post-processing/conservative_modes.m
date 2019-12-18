function [modes_PLS, E, C] = conservative_modes(modes, n_vars, n_atoms, D)
% This function creates the conservative modes and the modes reconstructed
% using the conservative modes and PLS prediction.
%
% If the elements on diagonal of modes_PLS' * modes are close to 1.0, the
% respective mode is a linear combination of the conservative modes.
%
% INPUT PARAMETERS --------------------------------------------------------
%
% - modes
%
%       a matrix with eigenvectors.
%
%       Eigvec-1     Eigvec-M
%     [                       ] weight 1
%     [                       ] weight 2
%     [                       ]   .
%     [         Modes         ]   .
%     [                       ]   .
%     [                       ]   .
%     [                       ] weight M
%
% - n_vars
%
%       number of variables.
%
% - n_atoms
%
%       number of atoms.
%
% - D
%
%       vector with variable scalings.
%
% OUTPUT PARAMETERS -------------------------------------------------------
%
% - modes_PLS
%
%       the reconstructed modes obtained from PLS regression. These modes 
%       are the linear combination of the conservative modes from matrix C.
%
% - E
%
%       matrix of mass fractions of each atom in each of the variables.
%
% - C
%
%       conservative modes.

C = zeros(n_vars, n_atoms);

C_m = 12.0107;
H_m = 1.00794;
O_m = 15.999;
N_m = 14.0067;

E = [%C H O N
    [ 0 1 0 0] % H2
    [0 0 1 0] % O2
    [0 0 1 0] % O
    [0 H_m/(O_m + H_m) O_m/(O_m + H_m) 0] % OH
    [0 (2*H_m)/(2*H_m + O_m) (O_m)/(2*H_m + O_m) 0] % H2O
    [0 1 0 0] % H
    [0 (H_m)/(H_m + 2*O_m) (2*O_m)/(H_m + 2*O_m) 0] % HO2
    [0 (2*H_m)/(2*H_m + 2*O_m) (2*O_m)/(2*H_m + 2*O_m) 0] % H2O2
    [C_m/(C_m + O_m) 0 O_m/(C_m + O_m) 0] % CO
    [C_m/(C_m + 2*O_m) 0 (2*O_m)/(C_m + 2*O_m) 0] % CO2
    [C_m/(H_m + C_m + O_m) H_m/(H_m + C_m + O_m) O_m/(H_m + C_m + O_m) 0] % HCO
    [0 0 0 1] % N2
    ];

% Construct conservative modes:
for i = 1:1:n_vars
    
   for j = 1:1:n_atoms
       
       C(i,j) = E(i,j) * D(i);

   end
   
end

modes_PLS = C * inv(C' * C) * C' * modes;

end