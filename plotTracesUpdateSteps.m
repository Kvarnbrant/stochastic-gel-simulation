function plotTracesUpdateSteps(nIterationsOuter, ...
                               nIterationsInner1, ...
                               nIterationsInner2, ...
                               mean_stats_inner)



step1 = nIterationsInner1/100;
step2 = nIterationsInner2/100;
step = step1+step2;
parameters = ["Mean degree","Mean angle","Mean length"];
for j = 1:3
    figure
    hold on
    indcs = nIterationsOuter/2:nIterationsOuter;
    for i = indcs
        x = ((i-1)*step+1):((i-1)*step+step1);
        plot(x,mean_stats_inner(x,j),'r')
        x = ((i-1)*step+1+step1):i*step;
        plot(x,mean_stats_inner(x,j),'g')
        if i == indcs(end)
            x = ((i-1)*step+1):((i-1)*step+step1);
            first = plot(x,mean_stats_inner(x,j),'r');
            x = ((i-1)*step+1+step1):i*step;
            second = plot(x,mean_stats_inner(x,j),'g');
        end
    end
    title(parameters(j));
    legend([first,second],'Update step 1','Update step 2')
    xticks(linspace((nIterationsOuter/2-1)*step+1,nIterationsOuter*step,3));
    labels=floor(linspace(nIterationsOuter/2,nIterationsOuter+1,3));
    xticklabels({num2str(labels(1)),num2str(labels(2)),num2str(labels(3))})
    xlabel('Outer iterations')
    ylabel(parameters(j))
end

end