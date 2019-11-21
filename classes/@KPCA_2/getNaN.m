function h=getNaN(obj,varargin)
%% Input Check:
if ~isa(obj,'cell') && ~isa(obj{1},'pcaData')
    error('Invalid input.');
end
nin = size(varargin,2);

%% Main:
cn=length(obj); %clusters number
h=zeros(cn,2);
temp={};
for j=1:cn
    temp{1}=obj{j}.a; temp{2}=obj{j}.ak;
    for i=1:2
        h(j,i)=any(any(isnan(temp{i})));
    end
end
[i,j]=find(h==1); h=[];
h=[i,j];

end