function prob = getProbability1(dist, ...
                               nConnections1,nConnections2, ...
                               newAngles,oldAngles, ...
                               d0,n0,a0, ...
                               c1,c2,c3, ...
                               connected)

% Acceptance probability for first update step

if (isempty(newAngles)) 
    newAngles = a0;
end

if (isempty(oldAngles))
    oldAngles = a0;
end

if connected
    delta_f1 = -c1*(dist-d0)^2;
    delta_f2 = c2*((nConnections1-1-n0)^2 + ...
                   (nConnections2-1-n0)^2 - ...
                   (nConnections1-n0)^2 - ...
                   (nConnections2-n0)^2);
    delta_f3 = c3*(sum((newAngles-a0).^2)-sum((oldAngles-a0).^2));

else
    delta_f1 = c1*(dist-d0)^2;
    delta_f2 = c2*((nConnections1+1-n0)^2 + ...
                   (nConnections2+1-n0)^2 - ...
                   (nConnections1-n0)^2 - ...
                   (nConnections2-n0)^2);
    delta_f3 = c3*(sum((newAngles-a0).^2)-sum((oldAngles-a0).^2));
end

delta_f = delta_f1 + delta_f2 + delta_f3;

hastings_ratio = exp(-delta_f);

prob = min(hastings_ratio,1);

end


