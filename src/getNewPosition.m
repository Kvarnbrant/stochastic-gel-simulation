function [newPos,q1,q2,ts] = getNewPosition(node,u,D,sigma)
% Proposes a new position newPos of the node 'node' along the line segment
% inside of the node region D with direction vector u. The new position is 
% sampled from a truncated normal distribution whose corresponding normal 
% distribution has standard deviation sigma. 

% Get the line segments in positive and negative u directions, starting
% from 'node'. 
[tPos,reachedLimitPos] = rayTraversal(node,u,D);
[tNeg,reachedLimitNeg] = rayTraversal(node,-u,D);

% If the line segment only consists of one point (i.e. if 'node' is a voxel
% corner and the voxels in positive and negative u direction are not in the
% node region), return the same node.
if (tNeg == 0) && (tPos == 0)
    newPos = node;
    q1 = 1;    % q1 = q2 in this case 
    q2 = 1;    % (the actual value doesn't matter, as long as it's non-zero.)
    ts = [0,0,0];
    return 
end

[nx, ny, nz] = size(D);
t = 0;
N1 = makedist('Normal','mu',0,'sigma',sigma);
truncN1 = truncate(N1,-tNeg,tPos);

% sample a point on the line segment, which is not allowed to be the same
% as newPos or one of the end points.
while t == 0 || t == -tNeg || t == tPos
    t = random(truncN1);   
end
ts=[-tNeg,tPos,t];


newPos = node + t*u; 
newPos = mod(newPos,[nx ny nz]); % periodic boundaries

q1 = pdf(truncN1,t);

% if the line segment inside of the node region is bounded (i.e. not
% infinite due to periodic boundaries), the same line segment is used when
% calculating q1 and q2
if (~reachedLimitNeg && ~reachedLimitPos)
    N2 = makedist('Normal','mu',t,'sigma',sigma);
    truncN2 = truncate(N2,-tNeg,tPos);
    q2 = pdf(truncN2,0);

% else, we need to find the line segment obtained when starting from newPos.
else
    t2Pos = rayTraversal(newPos,u,D);
    t2Neg = rayTraversal(newPos,-u,D);
    N2 = makedist('Normal','mu',0,'sigma',sigma);
    truncN2 = truncate(N2,-t2Neg,t2Pos);
    q2 = pdf(truncN2,-t);
end

end



