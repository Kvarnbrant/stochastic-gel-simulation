function [padding, realIndex] = getPadding( ...
                                                    intersections, ...
                                                    dist, ...
                                                    xMin, xMax, ...
                                                    yMin, yMax, ...
                                                    zMin, zMax ...
                                                    )

% Function taken from Nisslert 2004, but with minor adjustments

maxIndex = size(intersections,1);
padding = [];
realIndex = [];


% Note that the code assumes xMin = 0
index = find(intersections(:,1) >= xMax - dist); % changed from > to >=
realIndex = [realIndex;index];
padding = [padding; [intersections(index,1)-xMax intersections(index,2) intersections(index,3)]];

% x direction
index = find(intersections(:,1) < xMin + dist);
realIndex = [realIndex;index];
padding = [padding; [intersections(index,1)+xMax intersections(index,2) intersections(index,3)]];

% y direction 
tempVector = [intersections;padding];
index = find(tempVector(:,2) >= yMax - dist); % changed from > to >=
padding = [padding;[tempVector(index,1) tempVector(index,2)-yMax tempVector(index,3)]];

tempIndex = find(index > maxIndex);

if (~isempty(tempIndex))
    index(tempIndex) = realIndex(index(tempIndex)-maxIndex); 
end

realIndex = [realIndex;index];

tempVector = [intersections;padding];
index = find(tempVector(:,2) < yMin + dist); 
index2 = find(tempVector(index,2) >= 0); 
index = index(index2);  
padding = [padding;[tempVector(index,1) tempVector(index,2)+yMax tempVector(index,3)]];

tempIndex = find(index > maxIndex);
if (~isempty(tempIndex))
    index(tempIndex) = realIndex(index(tempIndex)-maxIndex);
end
realIndex = [realIndex;index];

% z direction
tempVector = [intersections;padding];
index = find(tempVector(:,3) >= zMax - dist);  % changed from > to >=
padding = [padding;[tempVector(index,1) tempVector(index,2) tempVector(index,3)-zMax]];

tempIndex = find(index > maxIndex);

if (~isempty(tempIndex))
    index(tempIndex) = realIndex(index(tempIndex)-maxIndex);
end

realIndex = [realIndex;index];

tempVector = [intersections;padding];
index = find(tempVector(:,3) < zMin + dist);  
index2 = find(tempVector(index,3) >= 0);  
index = index(index2);
padding = [padding;[tempVector(index,1) tempVector(index,2) tempVector(index,3)+zMax]];

tempIndex = find(index > maxIndex);

if (~isempty(tempIndex))
    index(tempIndex) = realIndex(index(tempIndex)-maxIndex);
end

realIndex = [realIndex;index];
