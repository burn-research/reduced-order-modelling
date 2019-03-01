% This function performs the rotation of the matrix of eigenvectors or
% eigenvalues according to the VARIMAX criterion of Kaiser. 
% The Varimax rotation (Kaiser, 1958) is the most popular orthogonal
% rotation. The axes in the new coordinate system are rotated so that the
% varimax criterion is maximized. The varimax criterion maximizes the sum
% of the variances of the squared loadings within each column of the
% loading matrix.

function [rotated_loadings] = rot(b, n_comp, n_obs, normalize)
if nargin < 4 || isempty(normalize)
    normalize = false;
end
if normalize
    hj = sqrt(diag(b*b'));
    Hj = repmat(hj, 1, size(b,2));
    loadings = b./Hj;
else
    loadings = b;
end

% Now we will set a maximum number of iteration, iter_max, and a relative
% tolerance, r_tol, for the judgment of convergence. The convergence will
% be reached when the amount of variance explained by each component
% loading will not vary from one iteration to another. If the maximum
% number of iteration is reached, iter_max, then the loop stops and an
% error messagen "No convergence, maximum number of iteration reached" will
% be displayed. 
   
iter_max = 2000;     % Maximum number of iteration (column pairs)
r_tol = 1e-16;      % Relative tolerance
iter = 1;           % Iteration inizialization
convergence = 0;    % Logical indicator of convergence
var = sum(loadings.^2);    % Variance explained

while (convergence == 0) && (iter < iter_max)
    for i = 1 : n_comp
        for j = i+1 : n_comp
            xj = loadings(:,i);
            yj = loadings(:,j);
            
            % This algorithm is based on the existence of an exact solution
            % for the VARIMAX criterion for the case n_comp=2 
            
            uj = xj.*xj - yj.*yj; 
            vj = 2*xj.*yj;
            A = sum(uj);
            B = sum(vj);
            C = uj'*uj - vj'*vj;
            D = 2*uj'*vj;
            num = D - 3*A*B/n_obs;
            den = C - (A^2 - B^2)/n_obs;
            tan4p = num/den;
                
            % Then we can evaluate the optimizing angle for rotation
            
            phi = atan2(num,den)/4;
            angle = phi*180/pi;
            if abs(phi) > 0.00001
                Xj = cos(phi)*xj + sin(phi)*yj;
                Yj = -sin(phi)*xj + cos(phi)*yj;
                bj1 = Xj;
                bj2 = Yj;
                b(:,i) = bj1;
                b(:,j) = bj2;
                loadings(:,i) = b(:,i);
                loadings(:,j) = b(:,j);
            end
        end
    end
    loadings = b;
        
    % We can test now the logical convergence flag "convergence".
     
    var_old = var;
    var = sum(loadings.^2);
    convergence = abs((var - var_old)/var) < r_tol;
    iter = iter + 1;
    
end

% If we performed Kaiser nromalization, the loadings must be renormalized
% to preserve the original variace explained. 

if normalize
    for i = 1 : n_comp
        rotated_loadings(:, i) = loadings(:, i) * Hj(i);
    end
else
    rotated_loadings = loadings;
end
% plot(loadings(:, 1), loadings(:, 2), 'bo');
end





