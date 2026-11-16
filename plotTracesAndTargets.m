function plotTracesAndTargets(nIterationsOuter,mean_stats_outer,mean_n,mean_angle,mean_dist, alt)

x = 0:nIterationsOuter; 
figure
subplot(2,2,1)
plot(x,mean_stats_outer(:,1),'b',x,mean_n*ones(nIterationsOuter+1),'r')
title('Trace of mean degree')
legend('Trace','Input target value')
xlabel('Number of outer iterations')
ylabel('Degree')

subplot(2,2,2)
plot(x,mean_stats_outer(:,2),'b',x,180/pi*mean_angle*ones(nIterationsOuter+1),'r')
title('Trace of mean angle between edges')
legend('Trace','Input target value')
xlabel('Number of outer iterations')
ylabel('Angle between edges')
subplot(2,2,3)

if alt % If using the alternative energy function, there is no target value of edge length
    plot(x,mean_stats_outer(:,3),'b')
    title('Trace of mean edge length')
    legend('Trace')
else
    plot(x,mean_stats_outer(:,3),'b',x,mean_dist*ones(nIterationsOuter+1),'r')
    title('Trace of mean edge length')
    legend('Trace','Input target value')
end

xlabel('Number of outer iterations')
ylabel('Edge length')

end
