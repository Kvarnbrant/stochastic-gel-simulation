function cumuPlot(mean_stats)

parameters = ["Mean degree","Mean angle","Mean edge length"];
[n, d] = size(mean_stats);
probs = [0.025 0.5 0.975];

for j = 1:d
    Q = NaN(n, length(probs));
    for t = 1:n
        Q(t, :) = quantile(mean_stats(1:t, j), probs);
    end
    
    figure
    hold on
    plot(1:n,Q(:,3),'--r','LineWidth',1.5)
    plot(1:n, Q(:,2),'b', 'LineWidth', 1.5)
    plot(1:n,Q(:,1),'--g','LineWidth',1.5)
    


    title('Cumulative quantile plot: '+[parameters(j)])
    legend('97.5% quantile','50% quantile','2.5% quantile')
    ylabel([parameters(j)])
    xlabel('Outer iterations')
end

end
