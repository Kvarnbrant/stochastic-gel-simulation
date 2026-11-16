function [X,Y,Z] = poissonPoints(nx,ny,nz,N)
% Input: size of the domain and number of points N.
% Output: Poisson distributed points in (0,nx)x(0,ny)x(0,nz) conditioned on
% N points.


X = nx*rand(N,1);
Y = ny*rand(N,1);
Z = nz*rand(N,1);


end
