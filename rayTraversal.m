function [t,reachedLimit] = rayTraversal(x0,u,D)

% Assumptions: 
% u isn't the zero vector
% x0 is inside the domain (both the box and the node region)
% x0 and u are row vectors
 
[nx,ny,nz]=size(D);
maxIt = nx + ny + nz;

u = u/norm(u);

i = floor(x0(1)) + 1;
j = floor(x0(2)) + 1;
k = floor(x0(3)) + 1;

step = sign(u);


nz_ind = find(u~=0);    % non-zero indices 
tDelta = zeros(1,3);   
tMax = inf(1,3);     % initialization

tDelta(nz_ind) = abs(1./u(nz_ind)); 

voxel_boundary = [i-1 + (step(1) > 0), j-1 + (step(2) > 0), k-1 + (step(3) > 0)];

tMax(nz_ind) = (voxel_boundary(nz_ind) - x0(nz_ind)) ./u(nz_ind);


reachedLimit = 0;
t = 0;  

for it = 1:maxIt
    tNext = min(tMax);
    mask = (tMax == tNext);

    t = tNext;
    i = i + mask(1)*step(1);
    j = j + mask(2)*step(2);
    k = k + mask(3)*step(3);

    tMax(mask) = tMax(mask) + tDelta(mask);

    % Periodic boundary wrapping
    if i > nx
        i = i - nx;
    elseif i < 1
        i = i + nx;
    end

    if j > ny
        j = j - ny;
    elseif j < 1
        j = j + ny;
    end

    if k > nz
        k = k - nz;
    elseif k < 1
        k = k + nz;
    end

    % if the limit of the node region is reached
    if ~D(i,j,k)    
        return
    end

end

reachedLimit = 1;

end


