function f1 = f1newEnergy(c1,n,nI,propTarget)


% If n==0, then f1 = c1(1)*propTarget(1) + c1(2)*propTarget(2) +
% c1(3)*propTarget(3)

if n == 0
    f1 = c1'*propTarget;
else
    f1 = c1(1)*abs(nI(1)/n-propTarget(1)) + c1(2)*abs(nI(2)/n-propTarget(2)) + c1(3)*abs(nI(3)/n-propTarget(3));
end

end

