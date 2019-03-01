function [opts, islocal] = setPcaOpts(this)

% PCA: Global or Local
t = this.pcaglobloc;
if isempty(t)
    t = input(['\nPCA: \n','[1] - Global \n','[2] - Local \n']);
end

islocal = false;
if t == 2
    islocal = true;
    opts.islocal = true;
end

% Correlation Domain
t = this.corrdomain;
if isempty(t)
    t = input('\nCorrelation Domain: [1] Spatial [2] Temporal: '); 
end

test = ismember(t,[1,2]);
if test
    opts.method = t;
else
    error('You had to choose <1> or <2>!');
end

% Constrained
t = this.con;
if isempty(t)
    t = input('\nConstrained? [y/n]    ');
end

opts.con = false;
if t == 'y' || t == 1 || t == true
    opts.con = true;
end

% Approximation Order
t = this.k;
if isempty(t)
    t = input('\Approximation Order: ');
end
opts.k = round(t); % App Ord has to be an integer

% Constraint function        
opts.mycon = @pcaData.allPositiveCon;

end


