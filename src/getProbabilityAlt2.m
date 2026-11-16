function prob = getProbabilityAlt2(I,deltaI,c1,totalEdges,propTarget,oldAngles,newAngles,propAngle,c3,q1,q2)

% No node degrees are changed in the second update step, hence delta_f2 = 0

f1Old = f1newEnergy(c1,totalEdges,I,propTarget);
f1New = f1newEnergy(c1,totalEdges,I+deltaI,propTarget);


delta_f3 =  c3*(sum((newAngles-propAngle).^2)-sum((oldAngles-propAngle).^2));

delta_f1 = f1New - f1Old;

delta_f = delta_f1 + delta_f3;

hastings_ratio = exp(-delta_f)*q2/q1;

prob = min(hastings_ratio,1);

end

