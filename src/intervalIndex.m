function distsIdx = intervalIndex(dists,l,r)

% Returns array containing 1's, 2's, 3's and 4's indicating the intervals
% that the lengths in dists belong to, where
% dists(j) in (0,l(1)) -> distsIdx(j) = 1
% dists(j) in [l(1),l(2)) -> distsIdx(j) = 2
% dists(j) in [l(2),r) -> distsIdx(j) = 3
% dists(j) in [r,inf] -> distsIdx(j) = 4

% Note that if dists contains edge lengths, distsIdx will never contain
% 4's, since no edge lengths can be longer than r. 

distsIdx = ones(size(dists)); 

idxInt2 = (dists>=l(1)) & (dists < l(2));
idxInt3 = (dists>=l(2)) & (dists < r);
idxInt4 = dists >= r;

distsIdx(idxInt2) = 2;
distsIdx(idxInt3) = 3;
distsIdx(idxInt4) = 4;


end

