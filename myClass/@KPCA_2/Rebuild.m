function Y = Rebuild(obj, varargin)
%% Explanation for input:
% varargin{1}=kri=1 --> rebuild the interpolated field
% varargin{2}=orig=1 --> rebuild the original field, even if it's CPCA

%%
nin = length(varargin);


%% Analizing the input:
% PCA or CPCA
l = obj.con;

% Default values
kri = false; 
orig = false; 
q = obj.k;


if nin > 0
    % User supplied value (Kriged data or not)
    kri = varargin{1};
    
    if nin > 1
        % User supplied value (alpha and not gamma, even if it is CPCA)
        orig = varargin{2};
        
        if nin > 2
            % User supplied approximation order
            q = varargin{3};
            % Check if the user-supplied q is too high
            if q > obj.k
                q = obj.k;
            end
        end % nin > 2
    end % nin > 1
end

% Get these local variables
dec.l = l;  dec.kri = kri;   dec.orig = orig;


%% Operating according to the input:
ai = decision(dec, obj);


%% Check for ai
if isempty(ai)
    Y = [];
    return;
end


%% Rebuilding:
modes = obj.modes(:, 1:q);
ai = ai(1:q, :);
t = size(ai, 2); 
Y = modes * ai;  

m = obj.m;      d = obj.d;
if obj.cs == false
    m = 0; d = 1;
end

for i = 1 : t
    Y(:,i) = m + Y(:,i) .* d; % Adds the mean and unscales
end

end


function ai = decision(dec, obj)
% dec.l: cpod or not;
% dec.kri: kriging or not;
% dec.orig: classic despite having a cpod;

p = {};  
temp = obj.a; 
temp_kriged = obj.a_kriged;

if obj.con 
    temp = obj.gamma; 
    temp_kriged = obj.gamma_kriged; 
end

%% Decision Matrix:
%%%%%% dec.orig = false
p{1,1,1} = obj.a;   p{1,2,1} = obj.a_kriged;         %kriging false: cpod false, cpod true
p{2,1,1} = temp;    p{2,2,1} = temp_kriged;          %kriging true: cpod false, cpod true
%%%%%% dec.orig = true
p{1,1,2} = obj.a;   p{1,2,2} = obj.a_kriged;         %kriging false: pod always
p{2,1,2} = obj.a;   p{2,2,2} = obj.a_kriged;         %kriging true: pod always

%% 
l = dec.l+1; k = dec.kri+1; o = dec.orig+1;   
ai = p{l,k,o};

end




