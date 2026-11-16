function [closestCoordsVec,dists] = closestCoordsVec(p, pts, nx, ny, nz)

% Given a point p and an array of points pts, returns the coordinates of
% the periodic images of each point in pts closest to p, as well as the
% distances. 

domain = [nx ny nz];

dx = pts - p;

dx = dx - domain.*floor(dx./domain + 0.5);
closestCoordsVec = p + dx;
dists = vecnorm(dx,2,2);

end

