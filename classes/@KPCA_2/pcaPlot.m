function pcaPlot(obj, varargin)
%%
%   Function used to create plots  


%% Set font size:
    setFigOpts;

    
%% Check inputs:
ivar = [];
nin = length(varargin); 
k = obj.k;
if nin > 0
    k = varargin{1}; 
    if ~(k==round(k)); error('2nd input: must be an integer.'); end
    if nin > 1
        var = varargin{2};
        if isa(var, 'char')
            ivar = strfind(var, obj.vars);
        else
            ivar = var;
        end
    end
end


%% figure(mode)
figure(); hold on;
for i = 1 : k
    name = ['PC', num2str(i)];
    if isempty(ivar)
        plot(obj.modes(:,i), 'DisplayName',name);
        xlabel('Variables');
    elseif strcmp(obj.x_status, 'variable')
        xs = (ivar - 1) * obj.nq + 1;
        xf = ivar * obj.nq;
        plot(obj.xq, obj.modes(xs:xf, i), 'DisplayName',name);
        xlabel('x [cm]');
    elseif strcmp(obj.x_status, 'parameter')
        error('Impossible to get a PC for one specified variable when space is a parameter.');
    end
end
grid on; hold off;
ylabel('PCs');
legend('show');



%% figure(a)
% iClust=1;
% [r, c] = size(a2plot); [~, ck] = size(ak2plot);
% psize=c+ck;
% p1=obj{1}.parmap; p2=obj{1}.parmapk;
% p=[p1;p2]; ps=sort(p);
% if isLocalPca
%     iClust=input(['\nWhich cluster to plot (type "0" for all of them)? [Total: ',...
%         num2str(max(obj{1}.idx)),']   ']);
% end
% if iClust==0
%     a=[]; ak=[];
%     for i=1:r
%         a=[a,a2plot(i,:)];
%         ak=[ak,ak2plot(i,:)];
%     end
%     p1=[1:length(a)];
%     p2=[1:length(ak)];
% else
%     a=a2plot(iClust,:); ak=ak2plot(iClust,:);
%     for i=1:psize
%         x=ps(i);
%         a(i)=pcaData.sampleit(x,p,[a,ak]);
%     end
%     p1=ps;
% end
% 
% figure(); plot(p1,a,'b'); grid on; 
% hold on; plot(p2,ak,'r'); hold off;
% xlabel(' '); ylabel(['PCA score: ',num2str(k)]);
% legend('PCA-scores','Interpolated PCA-scores');


end




function i = strfind(var, vars)

j = false; i = 0;
while(~j)
    i = i + 1;
    j = strcmp(var, vars{i});
end

end








