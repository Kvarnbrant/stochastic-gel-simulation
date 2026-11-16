function [bw1,bw2,bw3]=densityPlot(mean_stats)

sz = size(mean_stats,1);
samples = mean_stats(round(sz/2):end,1);
parameters = ["Mean degree","Mean angle","Mean edge length"];

figure;
subplot(2,2,1)
histogram(samples,'Normalization','pdf')
hold on;
[f, xi,bw1] = ksdensity(samples); 
plot(xi, f, 'LineWidth', 2); 
title([parameters(1),'Parameter Posterior Density']);
xlabel(parameters(1));
ylabel('Density');
grid on;

samples = mean_stats(round(sz/2):end,2);

subplot(2,2,2)
histogram(samples,'Normalization','pdf')
hold on;
[f, xi,bw2] = ksdensity(samples); 
plot(xi, f, 'LineWidth', 2);
title([parameters(2),'Parameter Posterior Density']);
xlabel(parameters(2));
ylabel('Density');
grid on;

samples = mean_stats(round(sz/2):end,3);

subplot(2,2,3)
histogram(samples,'Normalization','pdf')
hold on;
[f, xi,bw3] = ksdensity(samples); 
plot(xi, f, 'LineWidth', 2); 
title([parameters(3),'Parameter Posterior Density']);
xlabel(parameters(3));
ylabel('Density');
grid on;

end
