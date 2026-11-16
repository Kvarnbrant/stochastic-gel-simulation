function deltaI = deltaNodeI(newDists, oldDists, l)

% newDists and oldDists are the new and the old lengths, repsectively,
% of the edges that have changed. 

oldNodeI = zeros(3,1);
oldNodeI(1) = nnz(oldDists < l(1));
oldNodeI(2) = nnz((oldDists >= l(1)) & (oldDists < l(2)));
oldNodeI(3) = nnz(oldDists >= l(2));

newNodeI = zeros(3,1);
newNodeI(1) = nnz(newDists < l(1));
newNodeI(2) = nnz((newDists >= l(1)) & (newDists < l(2)));
newNodeI(3) = nnz(newDists >= l(2));


deltaI = newNodeI - oldNodeI;

end


