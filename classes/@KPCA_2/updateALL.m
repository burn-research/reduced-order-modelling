function obj = updateALL(obj)


% Check for str
if isempty(obj.str)
    fprintf('No update was queued.\n');
    return;
elseif ~isa(obj.str, 'char')
    error('Input is not a STRING.');
end


% Update PCA modes
if strcmp(obj.str,'k')
    update_modes(obj);
end


% Update PCA scores
if strcmp(obj.str,'k')
    update_a(obj);
end


% Update CPCA scores
if strcmp(obj.str,'k') || strcmp(obj.str,'con') || strcmp(obj.str,'mycon')
    update_gamma(obj);
end


% Update LPCA
C = {'clustering_dim'; 'num_clusts'; 'xp'};
if any( strcmp(obj.str, C) )
    if ~strcmp(obj.str,'xp')
        update_lpca(obj);
    end
    update_lpca_prop(obj);
end


% Update Kriging
C = {'k'; 'xp_kriged'; 'con'; 'trendFun'; 'corrFun'};
if any( strcmp(obj.str, C) ) 
    update_Kriging(obj);
end


% Update rebuilt and sorted data
% if ~isempty(obj.str)
%     obj = updateRebuiltAndSortedData(obj);
% end



% Clear the string
obj.str = [];


% New line
fprintf('\n');
end



