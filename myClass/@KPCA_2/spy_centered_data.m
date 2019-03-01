function spy_centered_data(obj, varargin)

n_args = length(varargin);
points = obj.training_points;
if n_args > 0 
    points = varargin{1};
end
if n_args > 1
    if iscell(varargin{2})
        var_list = varargin{2};
    elseif ischar(varargin{2})
        var_list = varargin(2);
    else
        var_list = obj.variable_names(varargin{2});
    end
else
    var_list = obj.variable_names;
end

n_vars = length(var_list);
y = zeros(n_vars,1);
n_points = size(points,1);
l = {};
Markers = { 'o' , '+' , '*' , 'x' };
LS = {'-', '--', '-.', ':'}; 
mk = 1;
figure();
for i = 1 : n_vars
    char_var = var_list{i};
    y = vectorize(obj.get_variable(points, char_var, 'centered', obj.is_mesh_variable));
    if n_points == 1 && size(obj.mesh,2) == 1
        plot(obj.mesh, y, 'LineStyle',LS{mk});%, 'Marker',Markers{mk});
    else
        plot(y, 'LineStyle',LS{mk});%, 'Marker',Markers{mk});
    end
    hold on;
    l{end+1} = char_var;
    % Marker
    mk = mk + 1;
    if mk > length(LS)
        mk = 1;
    end
end
legend(l);
grid on;
if n_points == 1 && size(obj.mesh,2) == 1
    xlabel('$ \mathbf{x \ [cm]} $','Interpreter','LaTex');
else
    xlabel('$ \mathbf{ } $','Interpreter','LaTex');
end
ylabel('$ \mathbf{Value \ [-]} $','Interpreter','LaTex');
tit = ['\textbf{Centered-scaled data; sc=',num2str(obj.scaling_criterion),'}'];
title(tit,'Interpreter','LaTex');

end




