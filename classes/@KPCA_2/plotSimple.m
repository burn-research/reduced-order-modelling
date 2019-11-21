function plotSimple(obj, varargin)

nin = length(varargin);
i = [];
if nin > 0
    i = varargin{1};
end
input_var = [];
if nin > 1
    input_var = varargin{2};
end


x_status = strcmp(obj.x_status,'variable');

if x_status
    prompt = '\nEnter case (column) you want to plot:  ';
    if isempty(i); 
        i = input(prompt); 
    end
    y_rec = obj.Y_sorted(:,i);     y = obj.Y_orig(:,i);
    if ~isempty(input_var)
        y_rec = getThisVar(y_rec, obj.map_rows, input_var);     
        y = getThisVar(y, obj.map_rows, input_var);
    end
    x_label = 'y';
else
    prompt = '\nEnter variable (row) you want to plot:  ';
    if isempty(i); i = input(prompt); end;
    y_rec = obj.Y_sorted(i,:);     y = obj.Y_orig(i,:);
    x_label = 'q';
end

subplot(1,2,1);
plot(y_rec); hold on; plot(y); grid on; hold off;
xlabel(x_label);
title('Compare');

subplot(1,2,2);
plot( abs(y_rec - y) ./ mean(abs(eps + y)) ); grid on;
xlabel(x_label);
title('Rel. diff.');

end


