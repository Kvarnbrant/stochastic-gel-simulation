function angles = getAngles(p1,neighborCoords)
% Returns a vector containing all the angles between edges of p1.

nNbs = size(neighborCoords,1);
nAngles = nNbs*(nNbs-1)*0.5;
angles = zeros(1,nAngles);
vecs = neighborCoords - p1;
norms = vecnorm(vecs,2,2); % 2-norm, dim = 2, i.e. norm per row
ind = 1;

for i =1:nNbs
    for j = i+1:nNbs
        cos_ang = dot(vecs(i,:),vecs(j,:)) / (norms(i) * norms(j));
        cos_ang = max(min(cos_ang,1),-1); % prevents acos error due to earlier rounding errors 
        angles(ind) = acos(cos_ang);
        ind = ind + 1;
    end
end

end

