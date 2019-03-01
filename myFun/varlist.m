function varlist(varargin)
   fprintf('Number of arguments: %d\n',nargin)
   celldisp(varargin)
   whos varargin
end