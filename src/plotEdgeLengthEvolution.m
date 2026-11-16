function plotEdgeLengthEvolution(histoNum,histoLengths,histoIdx)

figure;
hold on;
grid on;

% Colours progress from light blue to dark blue
lightColor = [0.70, 0.85, 1.00];
darkColor  = [0.00, 0.15, 0.45];

colors = [linspace(lightColor(1),darkColor(1),histoNum)', ...
          linspace(lightColor(2),darkColor(2),histoNum)', ...
          linspace(lightColor(3),darkColor(3),histoNum)'];

for i = 1:histoNum
    [f, xi] = ksdensity(histoLengths{i});
    plot(xi, f, 'LineWidth', 2);%, 'Color', colors(i,:));
end

title('Evolution of edge length distribution during the simulation')
legendStrings = string(histoIdx) + " iter.";
legend(legendStrings)
xlabel('Edge length')
ylabel('Probability density')

end