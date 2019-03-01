function pv = pv_B4(modes, eigs, n_pc)
% if iscell(modes)
%     n = length(modes);
%     pv = cell(n,1);
%     for ii = 1 : n
%         pv{ii} = pv_B4(modes{ii}, eigs{ii}, n_pc);
%     end
%     return
% end
[n_var, n_modes] = size(modes);
for ii = 1 : n_modes
    modes(:,ii) = modes(:,ii) * eigs(ii);
end
pv = [];
for ii = 1 : n_var
    [~, loc] = max( abs(modes(ii,:)) );
    if loc <= n_pc
        pv = [pv; ii];
    end
end
end

