function angles = getUpdatedAngles(Neighbor,NeighborInd,conn,nodeEdges,nodeCoords,neighborsIndex,neighborsCoords)

% Returns the angles between edges centered at Node, when the edge between 
% Node and Neighbor is shifted. If there are no angles centered at Node
% when the edge between Node and Neighbor is shifted, that is if Node then 
% has 0 or 1 edge, then an empty list is returned. 

connectedNodesIndex = find(nodeEdges);
if conn 
    connectedNodesIndex = connectedNodesIndex(connectedNodesIndex~=Neighbor);
end

nConnectedNodes = length(connectedNodesIndex);

if (nConnectedNodes > 0)
    [~,Locb] = ismember(connectedNodesIndex,neighborsIndex);
    if conn
        connectedNodesCoords = neighborsCoords(Locb,:);
    else
        connectedNodesCoords = zeros(nConnectedNodes+1,3);
        connectedNodesCoords(1:end-1,:) = neighborsCoords(Locb,:);
        connectedNodesCoords(end,:) = neighborsCoords(NeighborInd,:);
    end
    angles = getAngles(nodeCoords,connectedNodesCoords);
else
    angles = []; 
end

end


