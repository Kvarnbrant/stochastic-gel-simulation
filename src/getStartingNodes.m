function Nodes = getStartingNodes(D,nNodes)
% Samples the starting coordinates of the nNodes nodes inside of the node
% region D. 

% Note that the subscripts, in e.g. x direction, are integers 1,...,nx
% whereas voxel i in x direction lies within the half open interval [i-1,i).

[nx,ny,nz] = size(D);

ind = find(D);
Nodes = [];
nodesToGo = nNodes;
while nodesToGo > 0
    voxels = datasample(ind,nodesToGo);  % column vector, with replacement
    [X,Y,Z] = ind2sub([nx,ny,nz],voxels); % column vectors
    newNodes = [X,Y,Z] - 1 + rand(nodesToGo,3);
    Nodes = [Nodes;newNodes];
    Nodes = unique(Nodes,'rows');
    nodesToGo = nNodes - size(Nodes,1);
end



