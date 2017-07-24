function [leftT2,rightT2,isAdjust,dt] = AdjustPropSpeed_backup(leftT,rightT,dt0,CritK,KI1,KI2)
global TipStatesInv nActTip AllEle_global ConnList_global;
global Mat  epsKI;
leftT2 = leftT;
rightT2 = rightT;

isAdjust = 0;
dt = dt0;
maxKv = zeros(nActTip,1);
for i = 1 : nActTip
    num = TipStatesInv(i);
    if AllEle_global(num,9) < -5
       % continue;
    end
    if AllEle_global(num,10) > 1.1 &&  ConnList_global(num,3)*ConnList_global(num,4) > -1e-5
        %NF
            maxKv(i) = max([(KI1(i)/Mat.K1NF)^2 + (KI2(i)/Mat.K2NF)^2,CritK(i)]);
%         end
    else
        maxKv(i) = CritK(i);
    end
end
maxK  = max(maxKv);
alpha = 0.6;
coeff = (maxK-1);
if abs(coeff) > epsKI
    isAdjust = 1;
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
            dt = dt0*1.5;
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