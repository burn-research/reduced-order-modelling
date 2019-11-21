function [varargout] = createKpcaObj(obj, varargin)

this = KPCA_2(obj.Yas);

this.training_points = obj.p;
this.prediction_points = obj.p_miss;

I = ismember(obj.p_orig, obj.p_miss, 'rows');
this.original_data = obj.Y_orig(:,I);

% Output
varargout{1} = this;

end

