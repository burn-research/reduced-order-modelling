function ikeep = principal_variables(pca, varargin)
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

method = 'B2';  % default method
x = [];

switch nargin
  case 1
  case {2,3}
    if isa(varargin{1},'char')
      method = varargin{1};
      if nargin==3
        x = varargin{2};
      end
    else
      x = varargin{1};
      if nargin==3
        method = varargin{2};
      end
    end
  otherwise
    error('invalid use of PCA::principal_variables');
end

switch( upper(method) )
  
  case 'B2'   % B2 Method of Jolliffe (1972)
    
    nvar   = pca.nvar;
    neta   = pca.neta;
    eigVec = pca.Q;  % eigenvectors
    
    % set indices for discarded variables by looking at eigenvectors
    % corresponding to the discarded eigenvalues
    idisc = zeros(1,nvar-neta);
    for i=1:nvar-neta

      j = nvar-i+1;
      [evsrt,isrt] = sort( abs(eigVec(:,j)), 'descend' );
      
      % find the largest weight in this eigenvector
      % that has not yet been identified.
      for j=1:nvar
        ivar = isrt(j);
        if all( idisc ~= ivar )
          idisc(i) = ivar;
          break;
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
    
    if isempty(x)
      error('You must supply the data vector x when using the M2 method');
    end
    
    eta = x2eta(pca,x);  % the PCs based on the full set of x.
    
    nvarTot= pca.nvar;
    neta   = pca.neta;
    
    idiscard = [];
    j=1;
    q = nvarTot;
    while q>neta
      
      [npts,nvar] = size(x);
      m2cut = 1e12;
      
      for i=1:nvar
      
        % look at a PCA obtained from a subset of x.
        xs = [x(:,1:i-1), x(:,i+1:nvar)];
        pca2 = PCA( xs, pca.scaling, neta );
        etaSub = x2eta( pca2, xs );
        
        cov = etaSub'*eta;  % covariance of the two sets of PCs
        
        [U,S,V] = svd(cov,'econ');  % svd of the covariance
        m2 = trace( eta'*eta + etaSub'*etaSub - 2*S );
        
        if m2<m2cut
          m2cut=m2;
          idisc = i;
        end
      end
      
      % discard the selected variable
      x = [ x(:,1:idisc-1), x(:,idisc+1:nvar) ];

      % determine the original index for this discarded variable.
      ii = setdiff( 1:nvarTot, idiscard );
      idisc = ii(idisc);
      idiscard(j) = idisc;
      fprintf('discarding variable: %i\n',idisc);
      
      j=j+1;
      q = q-1;
      
    end
    
    ikeep = sort( setdiff( 1:nvarTot, idiscard ), 'ascend' );

    
  otherwise
    error( strcat('Invalid method "',method,'" for identifying principle variables') );
end
end




