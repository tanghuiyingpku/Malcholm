function dtnew = adjustTimeStep(dt,dD,dP)
w = 0.2;
dpmax = 1;
ddmax = 1e-3;
dtmin = 1e5;
for i = 1 : length(dP);
    dti = dt*(1+w)*ddmax/(abs(dD(i))+abs(dD(i+length(dP)))+w*ddmax);
    if dti < dtmin
        dtmin = dti;
    end
    dti = dt*(1+w)*dpmax/(abs(dP(i))+w*dpmax);
    if dti < dtmin
        dtmin = dti;
    end
end
dtnew = dtmin;
end
