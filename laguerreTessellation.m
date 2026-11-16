function T = laguerreTessellation(nx, ny, nz, ...
    number_of_regions_x, number_of_regions_y, number_of_regions_z, ...
    X, Y, Z, R)

% Written by RISE

% Check if the number of particles seems sufficient to not violate the
% periodicity condition in the distance calculations.
number_of_spheres = numel(X);
a = ( number_of_spheres / (nx * ny * nz) )^(1/3);
number_of_spheres_shortest_dim = min([nx, ny, nz]) * a;
if number_of_spheres_shortest_dim < 5
    warning(['Approximate number of spheres in shortest dimension is ' ...
        num2str(number_of_spheres_shortest_dim) ...
        '. Consider increasing number of spheres.']);
end

% Tessellation.
T = zeros(nx, ny, nz);

% Adjusted radii to facilitate formulating the "Laguerre distance" as a
% Euclidean distance.
R2_adjusted = sqrt(max(R)^2 - R .* R);

% Loop over regions.
for current_x_region = 1:number_of_regions_x
    for current_y_region = 1:number_of_regions_y
        for current_z_region = 1:number_of_regions_z

            % Lower and upper bounds of current region.
            x_min = round( (current_x_region - 1) / number_of_regions_x * ...
                nx + 1 );
            x_max = round( current_x_region / number_of_regions_x * nx );
            y_min = round( (current_y_region - 1) / number_of_regions_y * ...
                ny + 1 );
            y_max = round( current_y_region / number_of_regions_y * ny );
            z_min = round( (current_z_region - 1) / number_of_regions_z * ...
                nz + 1 );
            z_max = round( current_z_region / number_of_regions_z * nz );

            % Vectors for current region.
            x_vector = x_min:x_max;
            y_vector = y_min:y_max;
            z_vector = z_min:z_max;

            % Centers for current region.
            x_c = 0.5 * (x_min + x_max);
            y_c = 0.5 * (y_min + y_max);
            z_c = 0.5 * (z_min + z_max);

            % Shift circle coordinates.
            X_shift = X;
            Y_shift = Y;
            Z_shift = Z;

            X_shift(X_shift >= x_c + nx/2) = ...
                X_shift(X_shift >= x_c + nx/2) - nx;
            X_shift(X_shift <= x_c - nx/2) = ...
                X_shift(X_shift <= x_c - nx/2) + nx;
            Y_shift(Y_shift >= y_c + ny/2) = ...
                Y_shift(Y_shift >= y_c + ny/2) - ny;
            Y_shift(Y_shift <= y_c - ny/2) = ...
                Y_shift(Y_shift <= y_c - ny/2) + ny;
            Z_shift(Z_shift >= z_c + nz/2) = ...
                Z_shift(Z_shift >= z_c + nz/2) - nz;
            Z_shift(Z_shift <= z_c - nz/2) = ...
                Z_shift(Z_shift <= z_c - nz/2) + nz;

            % Detect spheres in region.
            idx_inside =  find( (X_shift >= x_min - 0.5) & ...
                (X_shift <= x_max - 0.5) & ...
                (Y_shift >= y_min - 0.5) & ...
                (Y_shift <= y_max - 0.5) & ...
                (Z_shift >= z_min - 0.5) & ...
                (Z_shift <= z_max - 0.5) );

            % Search edges of region.
            XX = [];
            YY = [];
            ZZ = [];
            for dim = 1:3
                for edge = 1:2
                    switch dim
                        case 1
                            switch edge
                                case 1
                                    xx = x_min;
                                case 2
                                    xx = x_max;
                            end
                            yy = y_vector;
                            zz = z_vector;
                        case 2
                            switch edge
                                case 1
                                    yy = y_min;
                                case 2
                                    yy = y_max;
                            end
                            xx = x_vector;
                            zz = z_vector;
                        case 3
                            switch edge
                                case 1
                                    zz = z_min;
                                case 2
                                    zz = z_max;
                            end
                            xx = x_vector;
                            yy = y_vector;
                    end

                    [XX_, YY_, ZZ_] = ndgrid(xx, yy, zz);
                    XX_ = XX_(:);
                    YY_ = YY_(:);
                    ZZ_ = ZZ_(:);

                    XX = cat(1, XX, XX_);
                    YY = cat(1, YY, YY_);
                    ZZ = cat(1, ZZ, ZZ_);

                end
            end
            XX = XX - 0.5;
            YY = YY - 0.5;
            ZZ = ZZ - 0.5;

            [~, IDX] = pdist2([X_shift Y_shift Z_shift R2_adjusted], ...
                [XX YY ZZ zeros(numel(XX), 1)], 'euclidean', ...
                'smallest', 1);

            idx_close = unique(IDX);

            % Index of all spheres relevant to current region.
            idx_relevant = union(idx_inside, idx_close);

            % Entire region.
            [XX, YY, ZZ] = ndgrid(x_vector, y_vector, z_vector);
            XX = XX - 0.5;
            YY = YY - 0.5;
            ZZ = ZZ - 0.5;
            XX = XX(:);
            YY = YY(:);
            ZZ = ZZ(:);

            [~, IDX] = pdist2([X_shift(idx_relevant), ...
                Y_shift(idx_relevant), Z_shift(idx_relevant), ...
                R2_adjusted(idx_relevant)], ...
                [XX, YY, ZZ, zeros(numel(XX), 1)], 'euclidean', ...
                'smallest', 1);

            % Add result to tessellation array.
            IDX = idx_relevant(IDX);
            IDX = reshape(IDX, [x_max - x_min + 1, ...
                y_max - y_min + 1, z_max - z_min + 1]);
            T(x_vector, y_vector, z_vector) = IDX;
        end
    end
end

