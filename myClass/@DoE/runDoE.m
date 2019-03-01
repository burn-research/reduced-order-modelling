function runDoE(obj)

if ~obj.lhs 
    % AS
    [obj.p, obj.Yas, obj.p_miss, obj.Pot, obj.Infl_rel] ...
        = fakeDoE(obj.p_orig, obj.Y_orig, obj.p0, obj.nSamplesMax, obj.cdp_filter);    
else 
    % DISCRETE LHS
    [obj.p, obj.p_miss] = discrete_lhs(obj.p_orig, obj.p0, obj.nSamplesMax);
    obj.Yas = sampleit(obj.p, obj.p_orig, obj.Y_orig);
    % Relative influece
%     obj.Infl_rel = samplingInfluence(obj.p, obj.Yas);
end

end

