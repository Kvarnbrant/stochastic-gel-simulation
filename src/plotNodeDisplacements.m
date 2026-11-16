function plotNodeDisplacements(tValues,acc2)

accTs = abs(tValues(acc2,3));
rejTs = abs(tValues(~acc2,3));
limTs = quantile(rejTs,0.99);

binedges = linspace(0,limTs);

figure
histogram(accTs(accTs<limTs),'Normalization','pdf','BinEdges',binedges,'FaceColor','#009E73');
hold on
histogram(rejTs(rejTs<limTs),'Normalization','pdf','BinEdges',binedges,'FaceColor','#E69F00');
title('Distributions of proposed node displacements');
legend('Accepted node displacements','Rejected node displacements')
xlabel('Proposed node displacement');
ylabel('Probability density');

end