function [CritTheta,CritK] = FindCritDir_Anisotropy(KI1,KI2)
global nTip nActTip nAct TipStatesInv;
global ConnList_global AllEle_global;
CritTheta = zeros(nTip,1);
CritK  = zeros(nTip,1);
%Find Growth Direction and Critical K calculation
for ii = 1 : nActTip
    num = TipStatesInv(ii);
    if num < 0.1
        f = 1;
    end
    [type,~] = TipType(num,AllEle_global,ConnList_global);
    a = -90;
    b = 90;
    %angle_ele = asind(AllEle_global(num,5));
    [Kmax,angle]= FindMaxTheta(a,b,KI1(ii),KI2(ii));
    CritK(ii) = Kmax;
    if type > 1.9
        angle = angle + 180;
    end
    CritTheta(ii) = angle;
end


end
