function [D, largestSize, prop] = largestPeriodicComponent(D)

% Finds the largest connected component of D, periodic boundaries
% considered

[nx, ny, nz] = size(D);
total = nnz(D);

if total == 0
    error('D is empty.');
end


% Finds ordinary connected components (6 connectivity)
CC = bwconncomp(D, 6);


componentSize = cellfun(@numel, CC.PixelIdxList);


% Creates label matrix
L = labelmatrix(CC);



% Examining boundaries to find periodic connections between components

% X periodicity
A = L(1,:,:);
B = L(nx,:,:);

mask = (A ~= 0) & (B ~= 0) & (A ~= B);

sx = double(A(mask));
tx = double(B(mask));

% Y periodicity
A = L(:,1,:);
B = L(:,ny,:);

mask = (A ~= 0) & (B ~= 0) & (A ~= B);

sy = double(A(mask));
ty = double(B(mask));


% Z periodicity
A = L(:,:,1);
B = L(:,:,nz);

mask = (A ~= 0) & (B ~= 0) & (A ~= B);

sz = double(A(mask));
tz = double(B(mask));


% Combines all periodic connections
s = [sx; sy; sz];
t = [tx; ty; tz];

clear A B mask sx tx sy ty sz tz

% Create graph with ordinary connected components as nodes and components
% that are periodically connected get an edge between them

if isempty(s)

    % No periodic connections at all.
    bins = (1:CC.NumObjects)';

else

    G = graph(s, t, [], CC.NumObjects);

    % Checks which ordinary connected components are periodically connected
    [bins, ~] = conncomp(G);

end

% Calculates the size of each periodic component
periodicSize = accumarray( ...
    bins(:), ...
    componentSize(:), ...
    [], ...
    @sum);


% Finds the largest periodic component
[largestSize, largestBin] = max(periodicSize);

keepComponents = find(bins == largestBin);

prop = largestSize / total;


% Update D to be the largest periodic component
D(:) = false;

for k = 1:numel(keepComponents)
    c = keepComponents(k);
    D(CC.PixelIdxList{c}) = true;
end


end

