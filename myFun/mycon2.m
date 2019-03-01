function [c,ceq]=mycon2(gamma,obj,m)
% c(x) is the array of nonlinear inequality constraints at x. fmincon attempts to satisfy
% c(x) <= 0 for all entries of c.
% ceq(x) is the array of nonlinear equality constraints at x. fmincon attempts to satisfy
% ceq(x) = 0 for all entries of ceq.
% For example, x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
Nvar=55; n=obj.nx/Nvar;

% roundedMean=round(obj.meanY,10);
% m=obj.meanY; m(m<0)=0; %To avoid problems with c(x)

%% Non-linear inequalities:
yave=(m);
temp=(obj.modes*gamma).*obj.d;
y=yave+temp;        

c = -min(y(3*n:end)); 

%% Non-linear equalities:
y=obj.m+(obj.modes*gamma).*obj.d;

t=zeros(n,Nvar);
for i=3:Nvar
    t(:,i-2)=y(1+(i-1)*n:i*n);
end

ceq = []; % max(abs(sum(t,2)-1));

end

