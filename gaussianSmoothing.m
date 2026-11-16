function D = gaussianSmoothing(S, sigma)

% Written by RISE

[nx, ny, nz] = size(S);

filter_size = max(2 * ceil(2 * sigma) + 1, 2 * ceil(1/8 * max([nx, ny, nz])) + 1);

D = imgaussfilt3(double(S), sigma, 'filtersize', filter_size, 'padding', 'circular');
% D = D / max(D(:));

end
