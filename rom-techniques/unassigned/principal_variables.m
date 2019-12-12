function ikeep = principal_variables(X, neta, eigVec, method, varargin)
%% Description
% Extract principal variables from a PCA
%
%Examples:
% ikeep = principal_variables( pca );
% ikeep = principal_variables( pca, 'B4' );
% ikeep = principal_variables( pca, X, 'M2' );
%
% INPUTS:
%  pca     - the PCA object.
%
%  x       - [OPTIONAL] data arranged with observations in rows and
%            variables in columns.  Note that this is only required for the
%            "M2" method. 
%  neta    - Number of PVs to find.
%  method  - [OPTIONAL] the method for determining the principal variables.
%            The following methods are currently supported:
%
%           "B4" : selects principal variables based on the variables
%                  contained in the eigenvectors corresponding to the
%                  largest eigenvalues.
%
%           "B2" : selects pvs based on variables contained in the smallest
%                  eigenvalues.  These are discarded and the remaining
%                  variables are used as the principal variables.  This is
%                  the default method.
%
%           "M2" : At each iteration, each remaining variable is analyzed
%                  via PCA.  This is a very expensive method.
%
% OUTPUTS:
%  ikeep - a vector of indices of retained variables.
%
% example:
%  p = PCA( x, 'vast', 5 );
%  ikeep = principal_variables(pca);

%% Input
if nargin < 4
    method = 'B2'; % default method
end

%% Recursive use
% Local PCA
if iscell(X) && iscell(eigVec)
    if (length(X) ~= length(eigVec))
        error('Cell-arrays of data matrices and PCA modes must have the same length.');
    end
    n_clusters = length(X);
    ikeep = cell(n_clusters,1);
    for ii = 1 : n_clusters
        ikeep{ii} = principal_variables(X{ii}, neta, eigVec{ii}, method);
    end
    return
end

%% Main
% [WARNING] 
% Only the B2 method works. The M2 method should.
%
switch( upper(method) )
  case 'B2'   % B2 Method of Jolliffe (1972)
    nvar = size(X,2); % nvar = pca.nvar;
    % Set indices for discarded variables by looking at eigenvectors
    % corresponding to the discarded eigenvalues
    idisc = zeros(1, nvar - neta);
    for i = 1 : nvar - neta
      j = nvar - i + 1;
      [~, isrt] = sort( abs(eigVec(:,j)), 'descend' );
      % find the largest weight in this eigenvector
      % that has not yet been identified.
      for j = 1 : nvar
        ivar = isrt(j);
        if all( idisc ~= ivar )
          idisc(i) = ivar;
          break
        end
      end
    end
    ikeep = sort( setdiff( 1:nvar, idisc ), 'ascend' );
    
  case 'B4'   % B4 Forward method
    nvar   = pca.nvar;
    neta   = pca.neta;
    eigVec = pca.Q;  % eigenvectors
    % set indices for retained variables by looking at eigenvectors
    % corresponding to the retained eigenvalues
    ikeep = zeros(1,neta);
    for i=1:neta
      [evsrt,isrt] = sort( abs(eigVec(:,i)), 'descend' );
      % find the largest weight in this eigenvector
      % that has not yet been identified.
      for j=1:nvar
        ivar = isrt(j);
        if all( ikeep ~= ivar )
          ikeep(i) = ivar;
          break;
        end
      end
    end
    ikeep = sort(ikeep,'ascend');
    
  case 'M2'   % Note: this is EXPENSIVE
    if isempty(X)
      error('You must supply the data matrix X when using the M2 method');
    end
%     eta = x2eta(pca, X); % the PCs based on the full set X
    eta = X * eigVec;
%     nvarTot = pca.nvar;
    nvarTot = size(X, 2);
%     neta = pca.neta; % [deleted row: neta is an input]
    idiscard = [];
    j = 1;
    q = nvarTot;
    while q > neta
      [~, nvar] = size(X);
      m2cut = 1e12;
      for i = 1 : nvar
        % Look at a PCA obtained from a subset of X
        XS = [X(:,1:i-1), X(:,i+1:nvar)];
%         pca2 = PCA(XS, pca.scaling, neta);
%         etaSub = x2eta(pca2, XS);
        [~, etaSub, ~] = pca(XS, 'Centered', true, 'Algorithm', 'svd');
        cov = etaSub' * eta;  % Covariance of the two sets of PCs
        [~, S, ~] = svd(cov, 'econ');  % SVD of the covariance
        m2 = trace(eta'*eta + etaSub'*etaSub - 2*S);
        if m2 < m2cut
          m2cut = m2;
          idisc = i;
        end
      end
      % Discard the selected variable
      X = [ X(:,1:idisc-1), X(:,idisc+1:nvar) ];
      % Determine the original index for this discarded variable
      ii = setdiff( 1:nvarTot, idiscard );
      idisc = ii(idisc);
      idiscard(j) = idisc;
      fprintf('discarding variable: %i\n', idisc);
      j = j + 1;
      q = q - 1;
    end
    ikeep = sort(setdiff(1:nvarTot, idiscard), 'ascend');
    
  otherwise
    error( strcat('Invalid method "',method,'" for identifying principle variables') );
end
end




