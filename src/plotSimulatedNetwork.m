function plotSimulatedNetwork(Nodes,Edges,nx,ny,nz,maxDist)

nNodes = size(Nodes,1);
figure
clf
hold on
box off
color1 = [242 173 46]/255;
color2 = [214 103 26]/255;

for i = 1:nNodes
    for j = i+1:nNodes
        if Edges(i,j)
            dist = norm(Nodes(i,:)-Nodes(j,:));
            if dist > maxDist
                newCoords = closestCoordsVec(Nodes(i,:), ...
                                             Nodes(j,:), ...
                                             nx,ny,nz);
                plot3([Nodes(i,1),newCoords(1)], ...
                      [Nodes(i,2),newCoords(2)], ...
                      [Nodes(i,3),newCoords(3)], ...
                      'color',color1, ...
                      'LineWidth',0.9)
                xlim([0-maxDist nx+maxDist]) 
                ylim([0-maxDist ny+maxDist])
                zlim([0-maxDist nz+maxDist])
                newCoords = closestCoordsVec(Nodes(j,:), ...
                                             Nodes(i,:), ...
                                             nx,ny,nz);
                plot3([Nodes(j,1),newCoords(1)], ...
                      [Nodes(j,2),newCoords(2)], ...
                      [Nodes(j,3),newCoords(3)], ...
                      'color',color1, ...
                      'LineWidth',0.9)
                xlim([0-maxDist nx+maxDist]) 
                ylim([0-maxDist ny+maxDist])
                zlim([0-maxDist nz+maxDist])
            else
                plot3([Nodes(i,1),Nodes(j,1)], ...
                      [Nodes(i,2),Nodes(j,2)], ...
                      [Nodes(i,3),Nodes(j,3)], ...
                      'color',color1, ...
                      'LineWidth',0.9)
                xlim([0-maxDist nx+maxDist]) 
                ylim([0-maxDist ny+maxDist])
                zlim([0-maxDist nz+maxDist])
            end
        end
    end
end

scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),5,color2,'filled')
axis equal
xlim([-10 nx+10]) 
ylim([-10 ny+10])
zlim([-10 nz+10])

view(3)

view(45,30)
hold off



end


