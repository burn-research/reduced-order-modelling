function plotSamples(obj)

[~, l] = size(obj.p0);
if l > 3
    printf('\nInput space dimension is higher than 3, I cannot plot.\n');
    return;
end

switch l
    case 1
        figure(); plotFun1(obj.p0);
        if ~isempty(obj.p)
            figure(); plotFun3(obj.p); hold on; 
        end
        if ~isempty(obj.p_miss)
            plotFun3(obj.p_miss);
        end
        hold off;
    case 2
        figure(); plotFun2(obj.p0);
        if ~isempty(obj.p)
            figure(); plotFun3(obj.p); hold on; 
        end
        if ~isempty(obj.p_miss)
            plotFun3(obj.p_miss);
        end
        hold off;
    case 3
        figure(); plotFun3(obj.p0);
        if ~isempty(obj.p)
            figure(); plotFun3(obj.p); hold on; 
        end
        if ~isempty(obj.p_miss)
            plotFun3(obj.p_miss);
        end
        hold off;
end

end



function plotFun1(x)
plot(x, 'o');
grid on;
end


function plotFun2(x)
plot(x(:,1), x(:,2), 'o');
grid on;
end


function plotFun3(x)
plot3(x(:,1), x(:,2), x(:,3), 'o');
grid on;
end





