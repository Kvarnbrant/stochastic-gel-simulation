function prob = getProbability2(oldDists,newDists,oldAngles,newAngles,d0,a0,c1,c3,q1,q2)

% Since no edges are shifted in update step 2, delta_f2 is always zero. 

delta_f1 = c1*(sum((newDists-d0).^2)-sum((oldDists-d0).^2));


delta_f3 =  c3*(sum((newAngles-a0).^2)-sum((oldAngles-a0).^2));

delta_f = delta_f1 + delta_f3;

hastings_ratio = exp(-delta_f)*q2/q1;

prob = min(hastings_ratio,1);

end


