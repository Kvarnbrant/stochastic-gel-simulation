function S = skeleton(T)

% Written by RISE

% Pad array with periodic boundary conditions to facilitate computing
% derivatives on the edges.
T_padded = padarray(T, [1 1 1], 'circular', 'both');

% Derivative operators.
d_dx_operator = [-1 ; 1];
d_dy_operator = reshape(d_dx_operator, [1 2 1]);
d_dz_operator = reshape(d_dx_operator, [1 1 2]);

% Compute derivatives.
dTdx = convn(T_padded, d_dx_operator, 'same');
dTdy = convn(T_padded, d_dy_operator, 'same');
dTdz = convn(T_padded, d_dz_operator, 'same');

% Identify points where any derivative is larger than zero. These points
% constitute the skeleton.
gradT = abs(dTdx) + abs(dTdy) + abs(dTdz);
gradT = gradT(2:end-1, 2:end-1, 2:end-1);

S = gradT > 0;

end
