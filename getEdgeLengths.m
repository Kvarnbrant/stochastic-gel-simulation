function lengths = getEdgeLengths(Nodes,Edges,nx,ny,nz,type)

boxSize = [nx,ny,nz];

[i, j] = find(triu(Edges, 1)); % get unique edges

d = Nodes(i,:) - Nodes(j,:);


d = d - boxSize .* round(d ./ boxSize);

distances = sqrt(sum(d.^2, 2));

if type == "all"
    lengths = distances;
elseif type == "max"
    lengths = max(distances);
elseif type == "min"
    lengths = min(distances);
elseif type == "mean"
    lengths = mean(distances);
end

end

