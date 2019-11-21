function plot_reconstruction_error(obj, y1, y2, title_string)

figure();

subplot(1,2,1);
plot(100 * y1); grid on;
xlabel('Observation index');
ylabel('Error [%]');
title(title_string);

subplot(1,2,2);
plot(100 * y2); grid on;
xlabel('Variable index');
ylabel('Error [%]');
title(title_string);

end

