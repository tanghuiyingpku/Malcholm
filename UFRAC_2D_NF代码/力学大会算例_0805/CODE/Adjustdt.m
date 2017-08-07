function  [leftT2,rightT2,dt]  = Adjustdt(leftT,rightT,dt0,Kmax)
leftT2 = leftT;
rightT2 = rightT;
epsKI = 0.1;
coeff = (Kmax-1);
if abs(coeff) >= epsKI
    if coeff >= epsKI % KIC is large
        if leftT < -1e5 ||  rightT > 1e5
            dt = dt0/1.5;
            rightT2 = dt0;
        end
        if rightT < 1e5 && leftT > -1e5
            dt = 0.5*(leftT+ rightT);
        end
        if rightT < 1e5 && dt0 < rightT
            rightT2 = dt0;
        end
    end
    if coeff <= -epsKI
        if leftT < -1e5 ||  rightT > 1e5
            leftT2 = dt0;
            dt = dt0*1.4;
        end
        if rightT < 1e5 && leftT > -1e5
            dt = 0.5*(leftT+ rightT);
        end
        if leftT > -1e5 && dt0 > leftT
            leftT2 = dt0;
        end
    end
    disp('Change Change Time Step');
end
end