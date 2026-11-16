function plotEdgeLengthsUpdateStep1(edgeLengths,nEdges,acc1)

addDists       = edgeLengths(acc1 ==  1);
remDists       = edgeLengths(acc1 == -1);
rejectAddDists = edgeLengths(acc1 ==  2);
rejectRemDists = edgeLengths(acc1 == -2);



figure;
hold on;
[f, xi] = ksdensity(addDists);
plot(xi, f, 'LineWidth', 2); 
hold on;


if length(addDists) > sum(nEdges)*2
    [f, xi] = ksdensity(addDists(sum(nEdges):end));
    plot(xi, f, 'LineWidth', 2); 
    hold on;
else
    [f, xi] = ksdensity(addDists(round(length(addDists)/2):end)); 
    plot(xi, f, 'LineWidth', 2);
    hold on;
end

[f, xi] = ksdensity(remDists); 
plot(xi, f, 'LineWidth', 2);

hold on;
[f, xi] = ksdensity(rejectAddDists);
plot(xi, f, 'LineWidth', 2);
hold on;
[f, xi] = ksdensity(rejectRemDists); 
plot(xi, f, 'LineWidth', 2);
title('Distributions of edge lengths for edge-shifting proposals');
xlabel('Edge length');
ylabel('Probability density');
grid on;

histogram(edgeLengths,'Normalization','pdf')
legend('All added edges', ...
       'Added edges (after burn-in)', ...
       'Removed edges','Rejected adding proposal', ...
       'Rejected removal proposal', ...
       'All proposals')

end